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
    currentIbor = args[3]
    futureIbor = args[4]
    fullPrice = args[5]

    px = self.fullPriceFromDM(settlementDate,
                              nextCoupon,
                              currentIbor,
                              futureIbor,
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
                 freqType: FinFrequencyTypes,
                 accrualType: FinDayCountTypes,
                 faceAmount: float = 100.0):
        ''' Create FinFloatingRateNote object given its maturity date, its
        quoted margin, coupon frequency, accrual type. Face is the size of
        the position and par is the notional on which price is quoted. '''

        checkArgumentTypes(self.__init__, locals())

        self._issueDate = issueDate
        self._maturityDate = maturityDate
        self._quotedMargin = quotedMargin
        self._freqType = freqType
        self._accrualType = accrualType
        self._flowDates = []
        self._frequency = FinFrequency(freqType)
        self._faceAmount = faceAmount   # This is the position size
        self._par = 100.0   # This is how price is quoted
        self._redemption = 1.0 # This is amount paid at maturity TODO NOT USED

        self._flowDates = []
        self._flowAmounts = []

        self._settlementDate = FinDate(1, 1, 1900)
        self._accruedInterest = None
        self._accruedDays = 0.0

        self._calculateFlowDates()

###############################################################################

    def _calculateFlowDates(self):
        ''' Determine the bond cashflow payment dates. '''

        # This should only be called once from init 

        calendarType = FinCalendarTypes.NONE
        busDayRuleType = FinBusDayAdjustTypes.NONE
        dateGenRuleType = FinDateGenRuleTypes.BACKWARD

        self._flowDates = FinSchedule(self._issueDate,
                                      self._maturityDate,
                                      self._freqType,
                                      calendarType,
                                      busDayRuleType,
                                      dateGenRuleType)._generate()

###############################################################################

    def fullPriceFromDM(self,
                        settlementDate: FinDate,
                        nextCoupon: float,  # The total reset coupon on NCD
                        currentIbor: float,  # Ibor discount to NCD
                        futureIbor: float,  # Future constant Ibor rates
                        dm: float):   # Discount margin
        ''' Calculate the full price of the bond from its discount margin (DM)
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. '''

        self.calcAccruedInterest(settlementDate, nextCoupon)

        dayCounter = FinDayCount(self._accrualType)

        q = self._quotedMargin
        numFlows = len(self._flowDates)
        
        # We discount using Libor over the period from settlement to the ncd
        (alpha, _, _) = dayCounter.yearFrac(settlementDate, self._ncd)
        df = 1.0 / (1.0 + alpha * (currentIbor + dm))

        # A full coupon is paid
        (alpha, _, _) = dayCounter.yearFrac(self._pcd, self._ncd)
        pv = nextCoupon * alpha * df
        
        # Now do all subsequent coupons that fall after the ncd
        for iFlow in range(1, numFlows):            

            if self._flowDates[iFlow] > self._ncd:

                pcd = self._flowDates[iFlow-1]
                ncd = self._flowDates[iFlow]
                (alpha, _, _) = dayCounter.yearFrac(pcd, ncd)

                df = df / (1.0 + alpha * (futureIbor + dm))
                c = futureIbor + q
                pv = pv + c * alpha * df

        pv += df
        pv = pv * self._par
        return pv

###############################################################################

    def principal(self,
                  settlementDate: FinDate,
                  nextCoupon: float,
                  currentIbor: float,
                  futureIbor: float,
                  dm: float):
        ''' Calculate the clean trade price of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. '''

        fullPrice = self.fullPriceFromDM(settlementDate,
                                         nextCoupon,
                                         currentIbor,
                                         futureIbor,
                                         dm)

        accrued = self._accruedInterest
        principal = fullPrice * self._faceAmount / self._par - accrued
        return principal

###############################################################################

    def dollarDuration(self,
                           settlementDate: FinDate,
                           nextCoupon: float,
                           currentIbor: float,
                           futureIbor: float,
                           dm: float):
        ''' Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. '''

        dy = 0.0001 # 1 basis point

        p0 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor + dy,
                                  futureIbor,
                                  dm)

        p2 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor - dy,
                                  futureIbor,
                                  dm)

        durn = (p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def dollarCreditDuration(self,
                             settlementDate: FinDate,
                             nextCoupon: float,
                             currentIbor: float,
                             futureIbor: float,
                             dm: float):
        ''' Calculate the risk or dP/dy of the bond by bumping. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calcAccruedInterest(settlementDate, nextCoupon)
        dy = 0.0001

        p0 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm + dy)

        p2 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm)

        durn = (p2 - p0) / dy
        return durn

###############################################################################

    def macauleyRateDuration(self,
                             settlementDate: FinDate,
                             nextCoupon: float,
                             currentIbor: float,
                             futureIbor: float,
                             dm: float):
        ''' Calculate the Macauley duration of the FRN on a settlement date
        given its yield to maturity. '''

        dd = self.dollarDuration(settlementDate,
                                 nextCoupon,
                                 currentIbor,
                                 futureIbor,
                                 dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm)

        md = dd * (1.0 + (nextCoupon + dm) / self._frequency) / fp
        return md

###############################################################################

    def modifiedRateDuration(self,
                             settlementDate: FinDate,
                             nextCoupon: float,
                             currentIbor: float,
                             futureIbor: float,
                             dm: float):
        ''' Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. '''

        dd = self.dollarDuration(settlementDate,
                                 nextCoupon,
                                 currentIbor,
                                 futureIbor,
                                 dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm)
        md = dd / fp
        return md

###############################################################################

    def modifiedCreditDuration(self,
                               settlementDate: FinDate,
                               nextCoupon: float,
                               currentIbor: float,
                               futureIbor: float,
                               dm: float):
        ''' Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. '''

        dd = self.dollarCreditDuration(settlementDate,
                                       nextCoupon,
                                       currentIbor,
                                       futureIbor,
                                       dm)

        fp = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm)
        md = dd / fp
        return md

###############################################################################

    def convexityFromDM(self,
                        settlementDate: FinDate,
                        nextCoupon: float,
                        currentIbor: float,
                        futureIbor: float,
                        dm: float):
        ''' Calculate the bond convexity from the discount margin (DM) using a
        standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. '''

        dy = 0.0001

        p0 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor - dy,
                                  futureIbor,
                                  dm)

        p1 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor,
                                  futureIbor,
                                  dm)

        p2 = self.fullPriceFromDM(settlementDate,
                                  nextCoupon,
                                  currentIbor + dy,
                                  futureIbor,
                                  dm)

        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def cleanPriceFromDM(self,
                         settlementDate: FinDate,
                         nextCoupon: float,
                         currentIbor: float,
                         futureIbor: float,
                         dm: float):
        ''' Calculate the bond clean price from the discount margin
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. '''

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        fullPrice = self.fullPriceFromDM(settlementDate,
                                         nextCoupon,
                                         currentIbor,
                                         futureIbor,
                                         dm)

        accrued = self._accruedInterest(settlementDate, nextCoupon)
        accrued = accrued * self._par / self._faceAmount

        cleanPrice = fullPrice - accrued
        return cleanPrice

###############################################################################

    def discountMargin(self,
                       settlementDate: FinDate,
                       nextCoupon: float,
                       currentIbor: float,
                       futureIbor: float,
                       cleanPrice: float):
        ''' Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. '''
        
        self.calcAccruedInterest(settlementDate, nextCoupon)

        # Needs to be adjusted to par notional
        accrued = self._accruedInterest * self._par / self._faceAmount

        fullPrice = cleanPrice + accrued

        argtuple = (self, settlementDate, nextCoupon, currentIbor,
                    futureIbor, fullPrice)

        dm = optimize.newton(_f,
                             x0=0.01,  # initial value of 10%
                             fprime=None,
                             args=argtuple,
                             tol=1e-12,
                             maxiter=50,
                             fprime2=None)

        return dm

###############################################################################

    def calcAccruedInterest(self,
                            settlementDate: FinDate,
                            nextCoupon: float):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Ex-dividend dates are 
        not handled. Contact me if you need this functionality. '''

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
                                          self._freqType)

        self._alpha = 1.0 - accFactor * self._frequency

        self._accruedInterest = accFactor * self._faceAmount * nextCoupon
        self._accruedDays = num
        return self._accruedInterest

###############################################################################

    def printFlows(self,
                   settlementDate: FinDate):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        self._calculateFlowDates()
        for dt in self._flowDates[1:-1]:
            print(dt)

        print(self._flowDates[-1])

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issueDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("QUOTED MARGIN (bp)", self._quotedMargin * 10000.0)
        s += labelToString("FREQUENCY", self._freqType)
        s += labelToString("ACCRUAL TYPE", self._accrualType)
        s += labelToString("FACE AMOUNT", self._faceAmount)
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
