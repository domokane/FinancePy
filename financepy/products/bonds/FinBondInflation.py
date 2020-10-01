###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# TODO: - ROUNDING CONVENTIONS FOR ACCRUED
# TODO: - CHECK OAS CALCULATION
# TODO:  - Check how first coupon on floating leg is sized on asset swaps. '''
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

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ..bonds.FinBond import FinBond, FinYTMCalcType

from scipy import optimize

###############################################################################


def _f(y, *args):
    ''' Function used to do root search in price to real yield calculation. '''
    bond = args[0]
    settlementDate = args[1]
    price = args[2]
    convention = args[3]
    px = bond.fullPriceFromRealYTM(settlementDate, y, convention)
    objFn = px - price
    return objFn

###############################################################################


def _g(oas, *args):
    ''' Function used to do root search in price to inflation yield. '''
    bond = args[0]
    settlementDate = args[1]
    price = args[2]
    baseCPI = args[3]
    referenceCPI = args[4]
    
    px = bond.fullPriceFromInflationYTM(settlementDate, 
                                        baseCPI,
                                        referenceCPI,
                                        )
    objFn = px - price
    return objFn

###############################################################################


class FinBondInflation(FinBond):
    ''' Class for inflation-linked bonds like TIPS and related analytics. These
    are bonds with coupon and principal adjusted by an index such as the CPI.
    '''

    def __init__(self,
                 issueDate: FinDate,
                 maturityDate: FinDate,
                 coupon: float,  # Annualised bond coupon before inflation
                 frequencyType: FinFrequencyTypes,
                 accrualType: FinDayCountTypes,
                 faceAmount: float = 100.0):
        ''' Create FinBond object by providing Maturity Date, Frequency,
        coupon and the accrual convention type. No issue date is required as
        we assume that the bond starts accruing before any settlement date and
        so that the first coupon is a full coupon. This may change. '''

        checkArgumentTypes(self.__init__, locals())

        if issueDate >= maturityDate:
            raise FinError("Issue Date must preceded maturity date.")

        self._issueDate = issueDate
        self._maturityDate = maturityDate
        self._coupon = coupon
        self._frequencyType = frequencyType
        self._accrualType = accrualType
        self._frequency = FinFrequency(frequencyType)
        self._faceAmount = faceAmount  # This is the bond holding size
        self._par = 100.0  # This is how price is quoted
        self._redemption = 1.0 # Amount paid at maturity

        self._flowDates = []
        self._flowAmounts = []

        self._settlementDate = FinDate(1, 1, 1900)
        self._accruedInterest = None
        self._accruedDays = 0.0
        self._alpha = 0.0

        # self._bond = FinBond(self._issueDate,
        #                      self._maturityDate,
        #                      self._coupon,
        #                      self._frequencyType,
        #                      self._accrualType,
        #                      self._faceAmount)
                   
        self._calculateFlowDates()
        self._calculateFlowAmounts()

# ###############################################################################

#     def _calculateFlowAmounts(self):
#         ''' Determine the bond cashflow payment amounts without principal and 
#         without any inflation adjustment. '''

#         self._flowAmounts = [0.0]

#         for dt in self._flowDates[1:]:
#            cpn = self._coupon / self._frequency
#            self._flowAmounts.append(cpn)
    
###############################################################################

    # def fullPriceFromRealYTM(self,
    #                          settlementDate: FinDate,
    #                          ytm: (float, np.ndarray),
    #                          convention: FinYTMCalcType=FinYTMCalcType.UK_DMO):
    #     ''' Calculate the full price of bond from its yield to maturity. This
    #     function is vectorised with respect to the yield input. It implements
    #     a number of standard conventions for calculating the YTM. '''

    #     bond = FinBond(self._issueDate, 
    #                    self._maturityDate, 
    #                    self._coupon,
    #                    self._frequencyType,
    #                    self._accrualType,
    #                    self._faceAmount)
                       
    #     price = bond.fullPriceFromYTM(settlementDate,
    #                                   ytm, 
    #                                   convention)

    #     return price

###############################################################################

    def fullPriceFromInflationYTM(self,
                                  settlementDate: FinDate,
                                  ytm: (float, np.ndarray),
                                  flatIndexRatio: float,
                                  indexRatio: float,
                                  convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the full price of bond from its inflation yield to
        maturity. This function is vectorised with respect to the yield input.
        It implements a number of standard conventions for calculating the YTM.
        '''

        bond = FinBond(self._issueDate, 
                       self._maturityDate, 
                       self._coupon * indexRatio,
                       self._frequencyType,
                       self._accrualType,
                       self._faceAmount)
                       
        fp = bond.fullPriceFromYTM(settlementDate,
                                   ytm, 
                                   convention)

        return fp * self._par

###############################################################################

    def principal(self,
                  settlementDate: FinDate,
                  ytm: float,
                  convention: FinYTMCalcType):
        ''' Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Libor rates. '''

        fullPrice = self.fullPriceFromYTM(settlementDate, ytm, convention)
        principal = fullPrice * self._faceAmount / self._par
        principal = principal - self._accruedInterest
        return principal

# ###############################################################################

#     def dollarDuration(self,
#                        settlementDate: FinDate,
#                        ytm: float,
#                        convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
#         ''' Calculate the risk or dP/dy of the bond by bumping. This is also
#         known as the DV01 in Bloomberg. '''

#         dy = 0.0001
#         p0 = self.fullPriceFromRealYTM(settlementDate, ytm - dy, convention)
#         p2 = self.fullPriceFromRealYTM(settlementDate, ytm + dy, convention)
#         durn = -(p2 - p0) / dy / 2.0
#         return durn

###############################################################################

#     def macauleyDuration(self,
#                          settlementDate: FinDate,
#                          ytm: float,
#                          convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
#         ''' Calculate the Macauley duration of the bond on a settlement date
#         given its yield to maturity. '''

#         dd = self.dollarDuration(settlementDate, ytm, convention)
#         fp = self.fullPriceFromRealYTM(settlementDate, ytm, convention)
#         md = dd * (1.0 + ytm / self._frequency) / fp
#         return md

# ###############################################################################

#     def modifiedDuration(self,
#                          settlementDate: FinDate,
#                          ytm: float,
#                          convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
#         ''' Calculate the modified duration of the bondon a settlement date
#         given its yield to maturity. '''

#         dd = self.dollarDuration(settlementDate, ytm, convention)
#         fp = self.fullPriceFromRealYTM(settlementDate, ytm, convention)
#         md = dd / fp
#         return md

# ###############################################################################

#     def convexityFromYTM(self,
#                          settlementDate: FinDate,
#                          ytm: float,
#                          convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
#         ''' Calculate the bond convexity from the yield to maturity. This
#         function is vectorised with respect to the yield input. '''

#         dy = 0.0001
#         p0 = self.fullPriceFromRealYTM(settlementDate, ytm - dy, convention)
#         p1 = self.fullPriceFromRealYTM(settlementDate, ytm, convention)
#         p2 = self.fullPriceFromRealYTM(settlementDate, ytm + dy, convention)
#         conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
#         return conv

###############################################################################

    def nominalYieldToMaturity(self, 
                               settlementDate: FinDate,
                               cleanPrice: float,
                               baseCPI: float,
                               refCPI: float,
                               convention: FinYTMCalcType=FinYTMCalcType.UK_DMO):
        
        indexRatio = refCPI / baseCPI
        cpn = self._coupon
        self._coupon = cpn * indexRatio
        self._redemption = indexRatio
        ytm = self.yieldToMaturity(settlementDate, cleanPrice, convention)
                                  
        return ytm

###############################################################################

    def cleanPriceFromNominalYTM(self,
                                 settlementDate: FinDate,
                                 ytm: float,
                                 convention: FinYTMCalcType=FinYTMCalcType.UK_DMO):
        ''' Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. '''

        
        fullPrice = self.fullPriceFromRealYTM(settlementDate, ytm, convention)
        accruedAmount = self._accruedInterest * self._par / self._faceAmount
        cleanPrice = fullPrice - accruedAmount
        return cleanPrice

###############################################################################

    def cleanPriceFromDiscountCurve(self,
                                    settlementDate: FinDate,
                                    discountCurve: FinDiscountCurve):
        ''' Calculate the clean bond value using some discount curve to
        present-value the bond's cashflows back to the curve anchor date and
        not to the settlement date. '''

        self.calcAccruedInterest(settlementDate)
        fullPrice = self.fullPriceFromDiscountCurve(settlementDate,
                                                    discountCurve)

        accrued = self._accruedInterest * self._par / self._faceAmount
        cleanPrice = fullPrice - accrued
        return cleanPrice

##############################################################################

    def fullPriceFromDiscountCurve(self,
                                   settlementDate: FinDate,
                                   discountCurve: FinDiscountCurve):
        ''' Calculate the bond price using a provided discount curve to PV the
        bond's cashflows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        '''

        if settlementDate < discountCurve._valuationDate:
            raise FinError("Bond settles before Discount curve date")

        if settlementDate > self._maturityDate:
            raise FinError("Bond settles after it matures.")

        px = 0.0
        df = 1.0
        dfSettle = discountCurve.df(settlementDate)

        for dt in self._flowDates[1:]:

            # coupons paid on the settlement date are included            
            if dt >= settlementDate:
                df = discountCurve.df(dt)
                flow = self._coupon / self._frequency
                pv = flow * df
                px += pv

        px += df
        px = px / dfSettle

        return px * self._par

###############################################################################

    def currentYield(self, cleanPrice):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._par / cleanPrice
        return y

###############################################################################

    def inflationYieldToMaturity(self,
                        settlementDate: FinDate,
                        cleanPrice: float,
                        convention: FinYTMCalcType=FinYTMCalcType.US_TREASURY):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''

        if type(cleanPrice) is float or type(cleanPrice) is np.float64:
            cleanPrices = np.array([cleanPrice])
        elif type(cleanPrice) is list or type(cleanPrice) is np.ndarray:
            cleanPrices = np.array(cleanPrice)
        else:
            raise FinError("Unknown type for cleanPrice "
                           + str(type(cleanPrice)))

        self.calcAccruedInterest(settlementDate)
        accruedAmount = self._accruedInterest * self._par / self._faceAmount
        fullPrices = (cleanPrices + accruedAmount)
        ytms = []

        for fullPrice in fullPrices:

            argtuple = (self, settlementDate, fullPrice, convention)

            ytm = optimize.newton(_g,
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

    def calcAccruedInterest(self, settlementDate: FinDate):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. '''

        numFlows = len(self._flowDates)

        if numFlows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        for iFlow in range(1, numFlows):
            # coupons paid on a settlement date are paid 
            if self._flowDates[iFlow] >= settlementDate:
                self._pcd = self._flowDates[iFlow-1]
                self._ncd = self._flowDates[iFlow]
                break
        
        dc = FinDayCount(self._accrualType)

        (accFactor, num, _) = dc.yearFrac(self._pcd,
                                          settlementDate,
                                          self._ncd, 
                                          self._frequency)

        self._alpha = 1.0 - accFactor * self._frequency
        self._accruedInterest = accFactor * self._faceAmount * self._coupon
        self._accruedDays = num
        
        return self._accruedInterest

###############################################################################



###############################################################################

    def fullPriceFromOAS(self,
                         settlementDate: FinDate,
                         discountCurve: FinDiscountCurve,
                         oas: float):
        ''' Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. '''

        self.calcAccruedInterest(settlementDate)
        f = self._frequency
        c = self._coupon

        pv = 0.0
        for dt in self._flowDates[1:]:
            
            # coupons paid on a settlement date are included
            if dt >= settlementDate:
                t = (dt - settlementDate) / gDaysInYear
                df = discountCurve.df(dt)
                # determine the Libor implied zero rate
                r = f * (np.power(df, -1.0 / t / f) - 1.0)
                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas)/f, -t * f)
                pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted
        pv *= self._par
        return pv

###############################################################################

    def printFlows(self,
                   settlementDate: FinDate):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''

        flow = self._faceAmount * self._coupon / self._frequency

        for dt in self._flowDates[1:-1]:
            # coupons paid on a settlement date are included
            if dt >= settlementDate:
                print("%12s" % dt, " %12.2f " % flow)

        redemptionAmount = self._faceAmount + flow
        print("%12s" % self._flowDates[-1], " %12.2f " % redemptionAmount)

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issueDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._frequencyType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount, "")
        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)


###############################################################################
