###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# TODO: - ROUNDING CONVENTIONS FOR ACCRUED
# TODO: - CHECK OAS CALCULATION
# TODO:  - Check how first coupon on floating leg is sized on asset swaps. """
###############################################################################

# https://www.dmo.gov.uk/media/15004/convention_changes.pdf

###############################################################################
# Conventions:
#  GILTS - SEMI ANNUAL ACT/ACT
#  US TREASURIES
###############################################################################

###############################################################################
# NOTE THAT I ASSUME THAT IF YOU SETTLE A SWAP ON A COUPON PAYMENT DATE YOU
# GET THE COUPON AND THE ACCRUED INTEREST EQUALS THE COUPON.
###############################################################################

import numpy as np

from ...utils.Date import Date
from ...utils.FinError import FinError
from ...utils.Frequency import Frequency, FinFrequencyTypes
from ...utils.FinGlobalVariables import gDaysInYear, gSmall
from ...utils.DayCount import DayCount, FinDayCountTypes
from ...utils.Schedule import Schedule
from ...utils.Calendar import Calendar
from ...utils.Calendar import FinCalendarTypes
from ...utils.Calendar import FinBusDayAdjustTypes
from ...utils.Calendar import FinDateGenRuleTypes
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from scipy import optimize

# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

###############################################################################


from enum import Enum


class FinYTMCalcType(Enum):
    UK_DMO = 1,
    US_STREET = 2,
    US_TREASURY = 3

###############################################################################


def _f(y, *args):
    """ Function used to do root search in price to yield calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    convention = args[3]
    px = bond.fullPriceFromYTM(settlement_date, y, convention)
    objFn = px - price
    return objFn

###############################################################################


def _g(oas, *args):
    """ Function used to do root search in price to OAS calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    discount_curve = args[3]
    px = bond.fullPriceFromOAS(settlement_date, discount_curve, oas)
    objFn = px - price
    return objFn

###############################################################################


class FinBond(object):
    """ Class for fixed coupon bonds and performing related analytics. These
    are bullet bonds which means they have regular coupon payments of a known
    size that are paid on known dates plus a payment of par at maturity. """

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,
                 coupon: float,  # Annualised bond coupon
                 freq_type: FinFrequencyTypes,
                 accrual_type: FinDayCountTypes,
                 face_amount: float = 100.0):
        """ Create FinBond object by providing the issue date, maturity Date,
        coupon frequency, annualised coupon, the accrual convention type, face
        amount and the number of ex-dividend days. """

        checkArgumentTypes(self.__init__, locals())

        if issue_date >= maturity_date:
            raise FinError("Issue Date must preceded maturity date.")

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._accrual_type = accrual_type
        self._frequency = Frequency(freq_type)
        self._face_amount = face_amount  # This is the bond holding size
        self._par = 100.0  # This is how price is quoted and amount at maturity
        self._redemption = 1.0 # This is amount paid at maturity

        self._flow_dates = []
        self._flow_amounts = []

        self._accruedInterest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

        self._calculateFlowDates()
        self._calculateFlowAmounts()

###############################################################################

    def _calculateFlowDates(self):
        """ Determine the bond cashflow payment dates."""

        # This should only be called once from init 

        calendar_type = FinCalendarTypes.NONE
        busDayRuleType = FinBusDayAdjustTypes.NONE
        date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

        self._flow_dates = Schedule(self._issue_date,
                                   self._maturity_date,
                                   self._freq_type,
                                   calendar_type,
                                   busDayRuleType,
                                   date_gen_rule_type)._generate()

###############################################################################

    def _calculateFlowAmounts(self):
        """ Determine the bond cashflow payment amounts without principal """

        self._flow_amounts = [0.0]

        for _ in self._flow_dates[1:]:
           cpn = self._coupon / self._frequency
           self._flow_amounts.append(cpn)
    
###############################################################################

    def fullPriceFromYTM(self,
                         settlement_date: Date,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the full price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM. """

        if convention not in FinYTMCalcType:
            raise FinError("Yield convention unknown." + str(convention))

        self.calcAccruedInterest(settlement_date)

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        f = Frequency(self._freq_type)
        c = self._coupon
        v = 1.0 / (1.0 + ytm/f)

        # n is the number of flows after the next coupon         
        n = 0
        for dt in self._flow_dates:
            if dt > settlement_date:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No coupons left")
 
        if convention == FinYTMCalcType.UK_DMO:
            if n == 0:
                fp = (v**(self._alpha))*(self._redemption + c/f)
            else:
                term1 = (c/f)
                term2 = (c/f)*v
                term3 = (c/f)*v*v*(1.0-v**(n-1))/(1.0-v)
                term4 = self._redemption * (v**n)
                fp = (v**(self._alpha))*(term1 + term2 + term3 + term4)
        elif convention == FinYTMCalcType.US_TREASURY:
            if n == 0:
                fp = (v**(self._alpha))*(self._redemption + c/f)
            else:
                term1 = (c/f)
                term2 = (c/f)*v
                term3 = (c/f)*v*v*(1.0-v**(n-1))/(1.0-v)
                term4 = self._redemption * (v**n)
                vw = 1.0 / (1.0 + self._alpha * ytm/f)
                fp = (vw)*(term1 + term2 + term3 + term4)
        elif convention == FinYTMCalcType.US_STREET:
            vw = 1.0 / (1.0 + self._alpha * ytm/f)
            if n == 0:
                vw = 1.0 / (1.0 + self._alpha * ytm/f)
                fp = vw*(self._redemption + c/f)
            else:
                term1 = (c/f)
                term2 = (c/f)*v
                term3 = (c/f)*v*v*(1.0-v**(n-1)) / (1.0-v)
                term4 = self._redemption * (v**n)
                fp = (v**(self._alpha))*(term1 + term2 + term3 + term4)
        else:
            raise FinError("Unknown yield convention")

        return fp * self._par

###############################################################################

    def principal(self,
                  settlement_date: Date,
                  y: float,
                  convention: FinYTMCalcType):
        """ Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. """

        fullPrice = self.fullPriceFromYTM(settlement_date, y, convention)

        principal = fullPrice * self._face_amount / self._par
        principal = principal - self._accruedInterest
        return principal

###############################################################################

    def dollarDuration(self,
                       settlement_date: Date,
                       ytm: float,
                       convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. """

        dy = 0.0001 # 1 basis point
        p0 = self.fullPriceFromYTM(settlement_date, ytm - dy, convention)
        p2 = self.fullPriceFromYTM(settlement_date, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def macauleyDuration(self,
                         settlement_date: Date,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. """

        dd = self.dollarDuration(settlement_date, ytm, convention)
        fp = self.fullPriceFromYTM(settlement_date, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

###############################################################################

    def modifiedDuration(self,
                         settlement_date: Date,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. """

        dd = self.dollarDuration(settlement_date, ytm, convention)
        fp = self.fullPriceFromYTM(settlement_date, ytm, convention)
        md = dd / fp
        return md

###############################################################################

    def convexityFromYTM(self,
                         settlement_date: Date,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dy = 0.0001
        p0 = self.fullPriceFromYTM(settlement_date, ytm - dy, convention)
        p1 = self.fullPriceFromYTM(settlement_date, ytm, convention)
        p2 = self.fullPriceFromYTM(settlement_date, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def cleanPriceFromYTM(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        """ Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        fullPrice = self.fullPriceFromYTM(settlement_date, ytm, convention)
        accruedAmount = self._accruedInterest * self._par / self._face_amount
        cleanPrice = fullPrice - accruedAmount
        return cleanPrice

###############################################################################

    def cleanPriceFromDiscountCurve(self,
                                    settlement_date: Date,
                                    discount_curve: FinDiscountCurve):
        """ Calculate the clean bond value using some discount curve to
        present-value the bond's cashflows back to the curve anchor date and
        not to the settlement date. """

        self.calcAccruedInterest(settlement_date)
        fullPrice = self.fullPriceFromDiscountCurve(settlement_date,
                                                    discount_curve)

        accrued = self._accruedInterest * self._par / self._face_amount
        cleanPrice = fullPrice - accrued
        return cleanPrice

##############################################################################

    def fullPriceFromDiscountCurve(self,
                                   settlement_date: Date,
                                   discount_curve: FinDiscountCurve):
        """ Calculate the bond price using a provided discount curve to PV the
        bond's cashflows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        """

        if settlement_date < discount_curve._valuation_date:
            raise FinError("Bond settles before Discount curve date")

        if settlement_date > self._maturity_date:
            raise FinError("Bond settles after it matures.")

        px = 0.0
        df = 1.0
        dfSettle = discount_curve.df(settlement_date)

        for dt in self._flow_dates[1:]:

            # coupons paid on the settlement date are included            
            if dt >= settlement_date:
                df = discount_curve.df(dt)
                flow = self._coupon / self._frequency
                pv = flow * df
                px += pv

        px += df * self._redemption
        px = px / dfSettle

        return px * self._par

###############################################################################

    def currentYield(self, cleanPrice):
        """ Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)"""

        y = self._coupon * self._par / cleanPrice
        return y

###############################################################################

    def yieldToMaturity(self,
                        settlement_date: Date,
                        cleanPrice: float,
                        convention: FinYTMCalcType = FinYTMCalcType.US_TREASURY):
        """ Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. """

        if type(cleanPrice) is float or type(cleanPrice) is np.float64:
            cleanPrices = np.array([cleanPrice])
        elif type(cleanPrice) is list or type(cleanPrice) is np.ndarray:
            cleanPrices = np.array(cleanPrice)
        else:
            raise FinError("Unknown type for cleanPrice "
                           + str(type(cleanPrice)))

        self.calcAccruedInterest(settlement_date)
        accruedAmount = self._accruedInterest * self._par / self._face_amount
        fullPrices = (cleanPrices + accruedAmount)
        ytms = []

        for fullPrice in fullPrices:

            argtuple = (self, settlement_date, fullPrice, convention)

            ytm = optimize.newton(_f,
                                  x0=0.05,  # guess initial value of 10%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        else:
            return np.array(ytms)

###############################################################################

    def calcAccruedInterest(self,
                            settlement_date: Date,
                            numExDividendDays: int = 0,
                            calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. If the
        bond trades with ex-coupon dates then you need to supply the number of
        days before the coupon date the ex-coupon date is. You can specify the
        calendar to be used - NONE means only calendar days, WEEKEND is only 
        weekends or you can specify a country calendar for business days."""

        numFlows = len(self._flow_dates)

        if numFlows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        for iFlow in range(1, numFlows):
            # coupons paid on a settlement date are paid 
            if self._flow_dates[iFlow] >= settlement_date:
                self._pcd = self._flow_dates[iFlow-1]
                self._ncd = self._flow_dates[iFlow]
                break

        dc = DayCount(self._accrual_type)
        cal = Calendar(calendar_type)
        exDividend_date = cal.addBusinessDays(self._ncd, -numExDividendDays)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd, 
                                            self._freq_type)
        
        if settlement_date > exDividend_date:
            acc_factor = acc_factor - 1.0 / self._frequency

        self._alpha = 1.0 - acc_factor * self._frequency
        self._accruedInterest = acc_factor * self._face_amount * self._coupon
        self._accrued_days = num
        
        return self._accruedInterest

###############################################################################

    def assetSwapSpread(
            self,
            settlement_date: Date,
            cleanPrice: float,
            discount_curve: FinDiscountCurve,
            swapFloatDayCountConventionType=FinDayCountTypes.ACT_360,
            swapFloatFrequencyType=FinFrequencyTypes.SEMI_ANNUAL,
            swapFloatCalendarType=FinCalendarTypes.WEEKEND,
            swapFloatBusDayAdjustRuleType=FinBusDayAdjustTypes.FOLLOWING,
            swapFloatDateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        """ Calculate the par asset swap spread of the bond. The discount curve
        is a Ibor curve that is passed in. This function is vectorised with
        respect to the clean price. """

        cleanPrice = np.array(cleanPrice)
        self.calcAccruedInterest(settlement_date)
        accruedAmount = self._accruedInterest * self._par / self._face_amount
        bondPrice = cleanPrice + accruedAmount
        # Calculate the price of the bond discounted on the Ibor curve
        pvIbor = 0.0
        prevDate = self._pcd

        for dt in self._flow_dates[1:]:
            
            # coupons paid on a settlement date are included
            if dt >= settlement_date:
                df = discount_curve.df(dt)
                pvIbor += df * self._coupon / self._frequency

        pvIbor += df * self._redemption

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prevDate = self._pcd
        schedule = Schedule(settlement_date,
                            self._maturity_date,
                            swapFloatFrequencyType,
                            swapFloatCalendarType,
                            swapFloatBusDayAdjustRuleType,
                            swapFloatDateGenRuleType)

        dayCount = DayCount(swapFloatDayCountConventionType)

        prevDate = self._pcd
        pv01 = 0.0
        for dt in schedule._adjustedDates[1:]:
            df = discount_curve.df(dt)
            year_frac = dayCount.year_frac(prevDate, dt)[0]
            pv01 = pv01 + year_frac * df
            prevDate = dt

        asw = (pvIbor - bondPrice/self._par) / pv01
        return asw

###############################################################################

    def fullPriceFromOAS(self,
                         settlement_date: Date,
                         discount_curve: FinDiscountCurve,
                         oas: float):
        """ Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. """

        self.calcAccruedInterest(settlement_date)
        f = self._frequency
        c = self._coupon

        pv = 0.0
        for dt in self._flow_dates[1:]:
            
            # coupons paid on a settlement date are included
            if dt >= settlement_date:

                t = (dt - settlement_date) / gDaysInYear

                t = np.maximum(t, gSmall)

                df = discount_curve.df(dt)
                # determine the Ibor implied zero rate
                r = f * (np.power(df, -1.0 / t / f) - 1.0)
                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas)/f, -t * f)
                pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted * self._redemption
        pv *= self._par
        return pv

###############################################################################

    def optionAdjustedSpread(self,
                             settlement_date: Date,
                             cleanPrice: float,
                             discount_curve: FinDiscountCurve):
        """ Return OAS for bullet bond given settlement date, clean bond price
        and the discount relative to which the spread is to be computed. """

        if type(cleanPrice) is float or type(cleanPrice) is np.float64:
            cleanPrices = np.array([cleanPrice])
        elif type(cleanPrice) is list or type(cleanPrice) is np.ndarray:
            cleanPrices = np.array(cleanPrice)
        else:
            raise FinError("Unknown type for cleanPrice "
                           + str(type(cleanPrice)))

        self.calcAccruedInterest(settlement_date)

        accruedAmount = self._accruedInterest * self._par / self._face_amount
        fullPrices = cleanPrices + accruedAmount

        oass = []

        for fullPrice in fullPrices:

            argtuple = (self, settlement_date, fullPrice, discount_curve)

            oas = optimize.newton(_g,
                                  x0=0.01,  # initial value of 1%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            oass.append(oas)

        if len(oass) == 1:
            return oass[0]
        else:
            return np.array(oass)

###############################################################################

    def printFlows(self,
                   settlement_date: Date):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        flow = self._face_amount * self._coupon / self._frequency

        for dt in self._flow_dates[1:-1]:
            # coupons paid on a settlement date are included
            if dt >= settlement_date:
                print("%12s" % dt, " %12.2f " % flow)

        redemptionAmount = self._face_amount + flow
        print("%12s" % self._flow_dates[-1], " %12.2f " % redemptionAmount)

###############################################################################

    def fullPriceFromSurvivalCurve(self,
                                   settlement_date: Date,
                                   discount_curve: FinDiscountCurve,
                                   survivalCurve: FinDiscountCurve,
                                   recovery_rate: float):
        """ Calculate discounted present value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the 
        defaulting principal we discretise the time steps using the coupon
        payment times. A finer discretisation may handle the time value with
        more accuracy. I reduce any error by averaging period start and period
        end payment present values. """

        f = self._frequency
        c = self._coupon

        pv = 0.0        
        prevQ = 1.0
        prevDf = 1.0

        defaultingPrincipalPVPayStart = 0.0
        defaultingPrincipalPVPayEnd = 0.0

        for dt in self._flow_dates[1:]:
            
            # coupons paid on a settlement date are included
            if dt >= settlement_date:

                df = discount_curve.df(dt)
                q = survivalCurve.survProb(dt)

                # Add PV of coupon conditional on surviving to payment date  
                # Any default results in all subsequent coupons being lost
                # with zero recovery

                pv = pv + (c / f) * df * q
                dq = q - prevQ

                defaultingPrincipalPVPayStart += -dq * recovery_rate * prevDf
                defaultingPrincipalPVPayStart += -dq * recovery_rate * df

                # Add on PV of principal if default occurs in coupon period
                prevQ = q
                prevDf = df

        pv = pv + 0.50 * defaultingPrincipalPVPayStart
        pv = pv + 0.50 * defaultingPrincipalPVPayEnd
        pv = pv + df * q * self._redemption
        pv *= self._par
        return pv

###############################################################################

    def cleanPriceFromSurvivalCurve(self,
                                    settlement_date: Date,
                                    discount_curve: FinDiscountCurve,
                                    survivalCurve: FinDiscountCurve,
                                    recovery_rate: float):
        """ Calculate clean price value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. """

        self.calcAccruedInterest(settlement_date)

        fullPrice = self.fullPriceFromSurvivalCurve(settlement_date,
                                                    discount_curve,
                                                    survivalCurve,
                                                    recovery_rate)
        
        cleanPrice = fullPrice - self._accruedInterest
        return cleanPrice
    
###############################################################################

    def __repr__(self):
        
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issue_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("ACCRUAL TYPE", self._accrual_type)
        s += labelToString("FACE AMOUNT", self._face_amount, "")
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)


###############################################################################
