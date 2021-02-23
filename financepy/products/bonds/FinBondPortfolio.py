###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# UNDER CONSTRUCTION !!!!!!!!!!!!!!!!

import numpy as np

from ...finutils.FinDate import FinDate
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from .FinBond import FinYTMCalcType

###############################################################################

class FinBondPortfolio(object):
    ''' Class for valuing and risk-managing a portfolio of bonds. '''

    def __init__(self,
                 bonds: (list),
                 bondWeights: (list, np.ndarray)):
        ''' XXX '''

        checkArgumentTypes(self.__init__, locals())

        self.calculateFlowAmounts()

###############################################################################

    def _calculateFlowAmounts(self):
        ''' Determine the bond cashflow payment amounts without principal '''

        self._flowAmounts = [0.0]

        for _ in self._flowDates[1:]:
           cpn = self._coupon / self._frequency
           self._flowAmounts.append(cpn)
    
###############################################################################

    def dollarDuration(self,
                       settlementDate: FinDate,
                       ytm: float,
                       convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. '''

        dy = 0.0001 # 1 basis point
        p0 = self.fullPriceFromYTM(settlementDate, ytm - dy, convention)
        p2 = self.fullPriceFromYTM(settlementDate, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def macauleyDuration(self,
                         settlementDate: FinDate,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm, convention)
        fp = self.fullPriceFromYTM(settlementDate, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

###############################################################################

    def modifiedDuration(self,
                         settlementDate: FinDate,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate, ytm, convention)
        fp = self.fullPriceFromYTM(settlementDate, ytm, convention)
        md = dd / fp
        return md

###############################################################################

    def convexityFromYTM(self,
                         settlementDate: FinDate,
                         ytm: float,
                         convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. '''

        dy = 0.0001
        p0 = self.fullPriceFromYTM(settlementDate, ytm - dy, convention)
        p1 = self.fullPriceFromYTM(settlementDate, ytm, convention)
        p2 = self.fullPriceFromYTM(settlementDate, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def cleanPriceFromYTM(self,
                          settlementDate: FinDate,
                          ytm: float,
                          convention: FinYTMCalcType = FinYTMCalcType.UK_DMO):
        ''' Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. '''

        fullPrice = self.fullPriceFromYTM(settlementDate, ytm, convention)
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


##############################################################################

    def fullPriceFromDiscountCurve(self,
                                   settlementDate: FinDate,
                                   discountCurve: FinDiscountCurve):
        ''' Calculate the bond price using a provided discount curve to PV the
        bond's cashflows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        '''


###############################################################################

    def currentYield(self, cleanPrice):
        ''' Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)'''

        y = self._coupon * self._par / cleanPrice
        return y

###############################################################################

    def yieldToMaturity(self,
                        settlementDate: FinDate,
                        cleanPrice: float,
                        convention: FinYTMCalcType = FinYTMCalcType.US_TREASURY):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''


###############################################################################

    def calcAccruedInterest(self, 
                            settlementDate: FinDate,
                            numExDividendDays: int = 0, 
                            calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND):
 
        return self._accruedInterest

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

    def fullPriceFromSurvivalCurve(self,
                                   settlementDate: FinDate,
                                   discountCurve: FinDiscountCurve,
                                   survivalCurve: FinDiscountCurve,
                                   recoveryRate: float):
        ''' Calculate discounted present value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the 
        defaulting principal we discretise the time steps using the coupon
        payment times. A finer discretisation may handle the time value with
        more accuracy. I reduce any error by averaging period start and period
        end payment present values. ''' 

        f = self._frequency
        c = self._coupon

        pv = 0.0        
        prevQ = 1.0
        prevDf = 1.0

        defaultingPrincipalPVPayStart = 0.0
        defaultingPrincipalPVPayEnd = 0.0

        for dt in self._flowDates[1:]:
            
            # coupons paid on a settlement date are included
            if dt >= settlementDate:

                df = discountCurve.df(dt)
                q = survivalCurve.survProb(dt)

                # Add PV of coupon conditional on surviving to payment date  
                # Any default results in all subsequent coupons being lost
                # with zero recovery

                pv = pv + (c / f) * df * q
                dq = q - prevQ

                defaultingPrincipalPVPayStart += -dq * recoveryRate * prevDf
                defaultingPrincipalPVPayStart += -dq * recoveryRate * df

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
                                   settlementDate: FinDate,
                                   discountCurve: FinDiscountCurve,
                                   survivalCurve: FinDiscountCurve,
                                   recoveryRate: float):
        ''' Calculate clean price value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. ''' 

        self.calcAccruedInterest(settlementDate)

        fullPrice = self.fullPriceFromSurvivalCurve(settlementDate,
                                                    discountCurve,
                                                    survivalCurve,
                                                    recoveryRate)
        
        cleanPrice = fullPrice - self._accruedInterest
        return cleanPrice
    
###############################################################################

    def __repr__(self):
        
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issueDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount, "")
        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)


###############################################################################
