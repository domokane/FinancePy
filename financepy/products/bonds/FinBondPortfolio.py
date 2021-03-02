###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# UNDER CONSTRUCTION !!!!!!!!!!!!!!!!

import numpy as np

from ...utils.Date import Date
from ...utils.Calendar import FinCalendarTypes
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from .FinBond import FinYTMCalcType

###############################################################################

class FinBondPortfolio(object):
    """ Class for valuing and risk-managing a portfolio of bonds. """

    def __init__(self,
                 bonds: (list),
                 bondWeights: (list, np.ndarray)):
        """ XXX """

        checkArgumentTypes(self.__init__, locals())

        self.calculateFlowAmounts()

###############################################################################

    def _calculateFlowAmounts(self):
        """ Determine the bond cashflow payment amounts without principal """

        self._flow_amounts = [0.0]

        for _ in self._flow_dates[1:]:
           cpn = self._coupon / self._frequency
           self._flow_amounts.append(cpn)
    
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


##############################################################################

    def fullPriceFromDiscountCurve(self,
                                   settlement_date: Date,
                                   discount_curve: FinDiscountCurve):
        """ Calculate the bond price using a provided discount curve to PV the
        bond's cashflows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        """


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


###############################################################################

    def calcAccruedInterest(self,
                            settlement_date: Date,
                            numExDividendDays: int = 0,
                            calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND):
 
        return self._accruedInterest

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
