##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from scipy import optimize


from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinCalendar import FinDateGenRuleTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################
# TODO: Need to complete and verify the risk sensitivity calculations.
###############################################################################


def _f(dm, *args):
    ''' Function used to do solve root search in DM calculation '''

    self = args[0]
    settlementDate = args[1]
    nextCoupon = args[2]
    currentLibor = args[3]
    futureLibor = args[4]
    fullPrice = args[5]

    px = self.fullPriceFromDM(settlementDate,
                              nextCoupon,
                              currentLibor,
                              futureLibor,
                              dm)

    objFn = px - fullPrice
    return objFn

###############################################################################


class FinBondFRN(object):
    ''' Class for managing floating rate notes that pay a floating index plus a
    quoted margin.'''

    def __init__(self,
                 issueDate: FinDate,
                 maturityDate: FinDate,
                 quotedMargin: float,    # Fixed spread paid on top of index
                 frequencyType: FinFrequencyTypes,
                 accrualType: FinDayCountTypes,
                 faceAmount: float = 100.0):
        ''' Create FinFloatingRateNote object given its maturity date, its
        quoted margin, coupon frequency, accrual type. Face is the size of
        the position and par is the notional on which price is quoted. '''

        print("Warning: An issue date argument will be added to this soon.")

        checkArgumentTypes(self.__init__, locals())

        self._issueDate = issueDate
        self._maturityDate = maturityDate
        self._quotedMargin = quotedMargin
        self._frequencyType = frequencyType
        self._accrualType = accrualType
        self._flowDates = []
        self._frequency = FinFrequency(frequencyType)
        self._faceAmount = faceAmount   # This is the position size
        self._par = 100.0   # This is how price is quoted

        ''' I do not determine cashflow dates as I do not want to require
        users to supply the issue date and without that I do not know how
        far to go back in the cashflow date schedule. '''

        self._settlementDate = FinDate(1900, 1, 1)
        self._accruedInterest = 0.0
        self._accruedDays = 0.0

        self._calculateFlowDates()

###############################################################################

    def _calculateFlowDates(self):
        ''' Determine the bond cashflow payment dates. '''

        calendarType = FinCalendarTypes.NONE
        busDayRuleType = FinBusDayAdjustTypes.NONE
        dateGenRuleType = FinDateGenRuleTypes.BACKWARD

        self._flowDates = FinSchedule(self._issueDate,
                                      self._maturityDate,
                                      self._frequencyType,
                                      calendarType,
                                      busDayRuleType,
                                      dateGenRuleType)._generate()

###############################################################################

    def fullPriceFromDM(self,
                        settlementDate: FinDate,
                        resetLibor: float,  # Next Libor payment on NCD
                        currentLibor: float,  # Libor discount to NCD
                        futureLibor: float,  # Future constant Libor rates
                        dm: float):   # Discount margin
        ''' Calculate the full price of the bond from its discount margin (DM)
        using standard model based on assumptions about future Libor rates. The
        next Libor payment which has reset is entered, so to is the current
        Libor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Libor payments and the discount
        margin. '''

        self.calcAccruedInterest(settlementDate, resetLibor)

        dayCounter = FinDayCount(self._accrualType)

        (_, dc, _) = dayCounter.yearFrac(settlementDate, self._ncd)
        (_, dbc, _) = dayCounter.yearFrac(self._pcd, self._ncd)

        alpha = float(dc) / float(dbc)
        f = self._frequency
        q = self._quotedMargin
        y = futureLibor + dm
        c = futureLibor + q
        numFlows = len(self._flowDates)
        dfPeriod = 1.0/(1.0 + y/f)
        
        # Handle the next coupon date (ncd)
        nextCoupon = resetLibor + self._quotedMargin
        df = pow(1.0 + (currentLibor + dm)/f, -alpha)
        pv = (nextCoupon/f) * df
        
        # Now do all subsequent coupons that fall after the ncd
        for iFlow in range(1, numFlows):            
            if self._flowDates[iFlow] > self._ncd:
                df *= dfPeriod
                pv = pv + (c/f) * df

        pv += df
        pv = pv * self._par
        return pv

###############################################################################

    def principal(self,
                  settlementDate: FinDate,
                  resetLibor: float,
                  currentLibor: float,
                  futureLibor: float,
                  dm: float):
        ''' Calculate the clean trade price of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Libor rates. '''

        fullPrice = self.fullPriceFromDM(settlementDate,
                                         resetLibor,
                                         currentLibor,
                                         futureLibor,
                                         dm)

        accrued = self._accruedInterest
        principal = fullPrice * self._faceAmount / self._par - accrued
        return principal

###############################################################################

    def dollarRateDuration(self,
                           settlementDate: FinDate,
                           resetLibor: float,
                           currentLibor: float,
                           futureLibor: float,
                           dm: float):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calcAccruedInterest(settlementDate, resetLibor)
        dy = 0.0001

        p0 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor + dy,
                                  futureLibor,
                                  dm)

        p2 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)

        durn = (p2 - p0) / dy
        return durn

###############################################################################

    def dollarCreditDuration(self,
                             settlementDate: FinDate,
                             resetLibor: float,
                             currentLibor: float,
                             futureLibor: float,
                             dm: float):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calcAccruedInterest(settlementDate, resetLibor)
        dy = 0.0001

        p0 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm + dy)

        p2 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)

        durn = (p2 - p0) / dy
        return durn

###############################################################################

    def macauleyRateDuration(self,
                             settlementDate: FinDate,
                             resetLibor: float,
                             currentLibor: float,
                             futureLibor: float,
                             dm: float):
        ''' Calculate the Macauley duration of the FRN on a settlement date
        given its yield to maturity. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dd = self.dollarRateDuration(settlementDate,
                                     resetLibor,
                                     currentLibor,
                                     futureLibor,
                                     dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)

        md = dd * (1.0 + (resetLibor + dm) / self._frequency) / fp
        return md

###############################################################################

    def modifiedRateDuration(self,
                             settlementDate: FinDate,
                             resetLibor: float,
                             currentLibor: float,
                             futureLibor: float,
                             dm: float):
        ''' Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Libor rates. The
        next Libor payment which has reset is entered, so to is the current
        Libor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Libor payments and the discount
        margin. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dd = self.dollarRateDuration(settlementDate,
                                     resetLibor,
                                     currentLibor,
                                     futureLibor,
                                     dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)
        md = dd / fp
        return md

###############################################################################

    def modifiedCreditDuration(self,
                               settlementDate: FinDate,
                               resetLibor: float,
                               currentLibor: float,
                               futureLibor: float,
                               dm: float):
        ''' Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Libor rates. The
        next Libor payment which has reset is entered, so to is the current
        Libor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Libor payments and the discount
        margin. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        dd = self.dollarCreditDuration(settlementDate,
                                       resetLibor,
                                       currentLibor,
                                       futureLibor,
                                       dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)
        md = dd / fp
        return md

###############################################################################

    def convexityFromDM(self,
                        settlementDate: FinDate,
                        resetLibor: float,
                        currentLibor: float,
                        futureLibor: float,
                        dm: float):
        ''' Calculate the bond convexity from the discount margin (DM) using a
        standard model based on assumptions about future Libor rates. The
        next Libor payment which has reset is entered, so to is the current
        Libor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Libor payments and the discount
        margin. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calcAccruedInterest(settlementDate, resetLibor)

        dy = 0.0001

        p0 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor - dy,
                                  futureLibor,
                                  dm)

        p1 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor,
                                  futureLibor,
                                  dm)

        p2 = self.fullPriceFromDM(settlementDate,
                                  resetLibor,
                                  currentLibor + dy,
                                  futureLibor,
                                  dm)

        print("The Value output has not been confirmed")
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def cleanPriceFromDM(self,
                         settlementDate: FinDate,
                         resetLibor: float,
                         currentLibor: float,
                         futureLibor: float,
                         dm: float):
        ''' Calculate the bond clean price from the discount margin
        using standard model based on assumptions about future Libor rates. The
        next Libor payment which has reset is entered, so to is the current
        Libor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Libor payments and the discount
        margin. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        fullPrice = self.fullPriceFromDM(settlementDate,
                                         resetLibor,
                                         currentLibor,
                                         futureLibor,
                                         dm)

        accrued = self._accruedInterest(settlementDate, resetLibor)
        accrued = accrued * self._par / self._faceAmount

        cleanPrice = fullPrice - accrued
        return cleanPrice

###############################################################################

    # def fullPriceFromDiscountCurve(self,
    #                                settlementDate: FinDate,
    #                                indexCurve: FinDiscountCurve,
    #                                discountCurve: FinDiscountCurve):
    #     ''' Calculate the bond price using some discount curve to present-value
    #     the bond's cashflows. THIS IS NOT COMPLETE. '''

    #     print("WARNING: DO NOT USE THIS FUNCTION")

    #     self._calculateFlowDates(settlementDate)

    #     pv = 0.0

    #     for dt in self._flowDates:
    #         df = discountCurve.df(dt)
    #         flow = self._coupon / self._frequency
    #         pv = pv + flow * df

    #     pv = pv + df
    #     return pv * self._faceAmount

###############################################################################

    def discountMargin(self,
                       settlementDate: FinDate,
                       resetLibor: float,
                       currentLibor: float,
                       futureLibor: float,
                       cleanPrice: float):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''
        
        self.calcAccruedInterest(settlementDate, resetLibor)

        # Needs to be adjusted to par notional
        accrued = self._accruedInterest * self._par / self._faceAmount

        fullPrice = cleanPrice + accrued

        argtuple = (self, settlementDate, resetLibor, currentLibor,
                    futureLibor, fullPrice)

        dm = optimize.newton(_f,
                             x0=0.01,  # initial value of 10%
                             fprime=None,
                             args=argtuple,
                             tol=1e-8,
                             maxiter=50,
                             fprime2=None)

        return dm

###############################################################################

    def calcAccruedInterest(self,
                            settlementDate: FinDate,
                            resetLibor: float):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. '''

        numFlows = len(self._flowDates)

        if numFlows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = FinDayCount(self._accrualType)

        for i in range(1, numFlows):
            if self._flowDates[i] > settlementDate:
                self._pcd = self._flowDates[i-1]
                self._ncd = self._flowDates[i]
                break

        (accFactor, num, _) = dc.yearFrac(self._pcd,
                                          settlementDate,
                                          self._ncd,
                                          self._frequencyType)

        self._alpha = 1.0 - accFactor * self._frequency
        nextCoupon = resetLibor + self._quotedMargin

        self._accruedInterest = accFactor * self._faceAmount * nextCoupon
        self._accruedDays = num
        return self._accruedInterest

###############################################################################

    def printFlows(self,
                   settlementDate: FinDate):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        self._calculateFlowDates(settlementDate)
        for dt in self._flowDates[1:-1]:
            print(dt)

        print(self._flowDates[-1])

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issueDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("QUOTED MARGIN (bp)", self._quotedMargin * 10000.0)
        s += labelToString("FREQUENCY", self._frequencyType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount)
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
