##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes, FinFrequency
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinCalendar, FinBusDayAdjustTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalTypes import FinSwapTypes

##########################################################################


class FinFixedIborXCcySwap(object):
    ''' Class for managing a cross currency swap contract. This is a contract
    in which a fixed payment leg in one currency is exchanged for floating
    payments in a second currency. There is an exchange of par at maturity.
    The contract is entered into at zero initial cost and it lasts from a start
    date to a specified maturity date.
    
    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curves (one for each leg) which is separate
    from the curve from which the implied index rates are extracted. '''
    
    def __init__(self,
                 effectiveDate: FinDate,  # Date interest starts to accrue
                 terminationDateOrTenor: (FinDate, str),  # Date contract ends
                 fixedLegType: FinSwapTypes,
                 fixedCoupon: float,  # Fixed coupon (annualised)
                 fixedFreqType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 floatSpread: float = 0.0,
                 floatFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create an interest rate swap contract giving the contract start
        date, its maturity, fixed coupon, fixed leg frequency, fixed leg day
        count convention and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. '''

        checkArgumentTypes(self.__init__, locals())

        if type(terminationDateOrTenor) == FinDate:
            self._terminationDate = terminationDateOrTenor
        else:
            self._terminationDate = effectiveDate.addTenor(terminationDateOrTenor)

        calendar = FinCalendar(calendarType)
        self._maturityDate = calendar.adjust(self._terminationDate,
                                             busDayAdjustType)

        if effectiveDate > self._maturityDate:
            raise FinError("Start date after maturity date")

        self._effectiveDate = effectiveDate
        self._notional = notional

        self._fixedCoupon = fixedCoupon
        self._floatSpread = floatSpread

        self._fixedFrequencyType = fixedFreqType
        self._floatFrequencyType = floatFreqType

        self._fixedDayCountType = fixedDayCountType
        self._floatDayCountType = floatDayCountType

        self._fixedLegType = fixedLegType

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        # These are generated immediately as they are for the entire
        # life of the swap. Given a valuation date we can determine
        # which cash flows are in the future and value the swap
        self._generateFixedLegPaymentDates()
        self._generateFloatLegPaymentDates()

        self._adjustedMaturityDate = self._adjustedFixedDates[-1]

        # Need to know latest payment date for bootstrap - DO I NEED THIS ??!
        self._lastPaymentDate = self._maturityDate
        if self._adjustedFixedDates[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjustedFixedDates[-1]

        if self._adjustedFloatDates[-1] > self._lastPaymentDate:
            self._lastPaymentDate = self._adjustedFloatDates[-1]

        # NOT TO BE PRINTED
        self._floatYearFracs = []
        self._floatFlows = []
        self._floatFlowPVs = []
        self._floatDfs = []

        self._fixedYearFracs = []
        self._fixedFlows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []

        self._firstFixingRate = None
        self._valuationDate = None
        self._fixedStartIndex = None

        self._calcFixedLegFlows()

##########################################################################

    def value(self,
              valuationDate,
              discountCurve,
              indexCurve,
              firstFixingRate=None,
              principal=0.0):
        ''' Value the interest rate swap on a value date given a single Libor
        discount curve. '''

        fixedLegValue = self.fixedLegValue(valuationDate,
                                           discountCurve,
                                           principal)

        floatLegValue = self.floatLegValue(valuationDate,
                                           discountCurve,
                                           indexCurve,
                                           firstFixingRate,
                                           principal)

        value = fixedLegValue - floatLegValue

        if self._fixedLegType == FinSwapTypes.PAY:
            value = value * (-1.0)

        return value

##########################################################################

    def _generateFixedLegPaymentDates(self):
        ''' Generate the fixed leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFixedDates = FinSchedule(
            self._effectiveDate,
            self._terminationDate,
            self._fixedFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType)._generate()

##########################################################################

    def _generateFloatLegPaymentDates(self):
        ''' Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFloatDates = FinSchedule(
            self._effectiveDate,
            self._terminationDate,
            self._floatFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType)._generate()

##########################################################################

    def fixedDates(self):
        ''' return a vector of the fixed leg payment dates '''
        if self._adjustedFixedDates is None:
            raise FinError("Fixed dates have not been generated")

        return self._adjustedFixedDates[1:]

##########################################################################

    def floatDates(self):
        ''' return a vector of the fixed leg payment dates '''
        if self._adjustedFloatDates is None:
            raise FinError("Float dates have not been generated")

        return self._adjustedFloatDates[1:]

##########################################################################

    def pv01(self, valuationDate, discountCurve):
        ''' Calculate the value of 1 basis point coupon on the fixed leg. '''

        pv = self.fixedLegValue(valuationDate, discountCurve)
        pv01 = pv / self._fixedCoupon / self._notional
        return pv01

##########################################################################

    def swapRate(self, valuationDate, discountCurve):
        ''' Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. '''

        pv01 = self.pv01(valuationDate, discountCurve)

        if valuationDate < self._effectiveDate:
            df0 = discountCurve.df(self._effectiveDate)
        else:
            df0 = discountCurve.df(valuationDate)

        dfT = discountCurve.df(self._maturityDate)

        if abs(pv01) < gSmall:
            raise FinError("PV01 is zero. Cannot compute swap rate.")

        cpn = (df0 - dfT) / pv01
        return cpn

##########################################################################

    def fixedLegValue(self, valuationDate, discountCurve, principal=0.0):

        self._valuationDate = valuationDate

        self._fixedYearFracs = []
        self._fixedFlows = []
        self._fixedDfs = []
        self._fixedFlowPVs = []
        self._fixedTotalPV = []

        dayCounter = FinDayCount(self._fixedDayCountType)

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. '''
        startIndex = 0
        while self._adjustedFixedDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._effectiveDate:
            startIndex = 1

        self._fixedStartIndex = startIndex

        ''' Now PV fixed leg flows. '''
        self._dfValuationDate = discountCurve.df(valuationDate)

        pv = 0.0
        prevDt = self._adjustedFixedDates[startIndex - 1]
        df_discount = 1.0
        if len(self._adjustedFixedDates) == 1:
            return 0.0

        for nextDt in self._adjustedFixedDates[startIndex:]:
            alpha = dayCounter.yearFrac(prevDt, nextDt)[0]
            df_discount = discountCurve.df(nextDt) / self._dfValuationDate
            flow = self._fixedCoupon * alpha * self._notional
            flowPV = flow * df_discount
            pv += flowPV
            prevDt = nextDt

            self._fixedYearFracs.append(alpha)
            self._fixedFlows.append(flow)
            self._fixedDfs.append(df_discount)
            self._fixedFlowPVs.append(flow * df_discount)
            self._fixedTotalPV.append(pv)

        flow = principal * self._notional
        pv = pv + flow * df_discount
        self._fixedFlowPVs[-1] += flow * df_discount
        self._fixedFlows[-1] += flow
        self._fixedTotalPV[-1] = pv
        return pv

##########################################################################

    def _calcFixedLegFlows(self):

        self._fixedYearFracs = []
        self._fixedFlows = []

        dayCounter = FinDayCount(self._fixedDayCountType)

        ''' Now PV fixed leg flows. '''
        prevDt = self._adjustedFixedDates[0]

        for nextDt in self._adjustedFixedDates[1:]:
            alpha = dayCounter.yearFrac(prevDt, nextDt)[0]
            flow = self._fixedCoupon * alpha * self._notional
            prevDt = nextDt
            self._fixedYearFracs.append(alpha)
            self._fixedFlows.append(flow)

##########################################################################

    def cashSettledPV01(self,
                        valuationDate,
                        flatSwapRate,
                        frequencyType):
        ''' Calculate the forward value of an annuity of a forward starting
        swap using a single flat discount rate equal to the swap rate. This is
        used in the pricing of a cash-settled swaption in the FinIborSwaption
        class. This method does not affect the standard valuation methods.'''

        m = FinFrequency(frequencyType)

        if m == 0:
            raise FinError("Frequency cannot be zero.")

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. '''
        startIndex = 0
        while self._adjustedFixedDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._effectiveDate:
            startIndex = 1

        ''' Now PV fixed leg flows. '''
        flatPV01 = 0.0
        df = 1.0
        alpha = 1.0 / m

        for _ in self._adjustedFixedDates[startIndex:]:
            df = df / (1.0 + alpha * flatSwapRate)
            flatPV01 += df * alpha

        return flatPV01

##########################################################################

    def floatLegValue(self,
                      valuationDate,  # This should be the settlement date
                      discountCurve,
                      indexCurve,
                      firstFixingRate=None,
                      principal=0.0):
        ''' Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. The valuation date can
        be the today date. In this case the price of the floating leg will not
        be par (assuming we added on a principal repayment). This is only the
        case if we set the valuation date to be the swap's actual settlement
        date. '''

        self._valuationDate = valuationDate
        self._floatYearFracs = []
        self._floatFlows = []
        self._floatRates = []
        self._floatDfs = []
        self._floatFlowPVs = []
        self._floatTotalPV = []
        self._firstFixingRate = firstFixingRate

        basis = FinDayCount(self._floatDayCountType)

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the start date. '''
        startIndex = 0
        while self._adjustedFloatDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._effectiveDate:
            startIndex = 1

        self._floatStartIndex = startIndex

        # Forward price to settlement date (if valuation is settlement date)
        self._dfValuationDate = discountCurve.df(valuationDate)

        ''' The first floating payment is usually already fixed so is
        not implied by the index curve. '''
        prevDt = self._adjustedFloatDates[startIndex - 1]
        nextDt = self._adjustedFloatDates[startIndex]
        alpha = basis.yearFrac(prevDt, nextDt)[0]
        df1_index = indexCurve.df(self._effectiveDate)  # Cannot be pcd as has past
        df2_index = indexCurve.df(nextDt)

        floatRate = 0.0

        if self._firstFixingRate is None:
            fwdRate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwdRate + self._floatSpread) * alpha * self._notional
            floatRate = fwdRate
        else:
            flow = self._firstFixingRate * alpha * self._notional
            floatRate = self._firstFixingRate

        # All discounting is done forward to the valuation date
        df_discount = discountCurve.df(nextDt) / self._dfValuationDate

        pv = flow * df_discount

        self._floatYearFracs.append(alpha)
        self._floatFlows.append(flow)
        self._floatRates.append(floatRate)
        self._floatDfs.append(df_discount)
        self._floatFlowPVs.append(flow * df_discount)
        self._floatTotalPV.append(pv)

        prevDt = nextDt
        df1_index = indexCurve.df(prevDt)

        for nextDt in self._adjustedFloatDates[startIndex + 1:]:
            alpha = basis.yearFrac(prevDt, nextDt)[0]
            df2_index = indexCurve.df(nextDt)
            # The accrual factors cancel
            fwdRate = (df1_index / df2_index - 1.0) / alpha
            flow = (fwdRate + self._floatSpread) * alpha * self._notional

            # All discounting is done forward to the valuation date
            df_discount = discountCurve.df(nextDt) / self._dfValuationDate

            pv += flow * df_discount
            df1_index = df2_index
            prevDt = nextDt

            self._floatFlows.append(flow)
            self._floatYearFracs.append(alpha)
            self._floatRates.append(fwdRate)
            self._floatDfs.append(df_discount)
            self._floatFlowPVs.append(flow * df_discount)
            self._floatTotalPV.append(pv)

        flow = principal * self._notional
        pv = pv + flow * df_discount
        self._floatFlows[-1] += flow
        self._floatFlowPVs[-1] += flow * df_discount
        self._floatTotalPV[-1] = pv

        return pv

##########################################################################

    def printFixedLegPV(self):
        ''' Prints the fixed leg dates, accrual factors, discount factors,
        cash amounts, their present value and their cumulative PV using the
        last valuation performed. '''

        print("START DATE:", self._effectiveDate)
        print("MATURITY DATE:", self._maturityDate)
        print("COUPON (%):", self._fixedCoupon * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixedFrequencyType))
        print("FIXED LEG DAY COUNT:", str(self._fixedDayCountType))
        print("VALUATION DATE", self._valuationDate)

        if len(self._fixedFlows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_DATE     YEAR_FRAC        FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        if self._fixedStartIndex is None:
            raise FinError("Need to value swap before calling this function.")

        startIndex = self._fixedStartIndex

        # By definition the discount factor is 1.0 on the valuation date
        print("%15s %10s %12s %12.8f %12s %12s" %
              (self._valuationDate,
               "-",
               "-",
               1.0,
               "-",
               "-"))

        iFlow = 0
        for paymentDate in self._adjustedFixedDates[startIndex:]:
            print("%15s %10.7f %12.2f %12.8f %12.2f %12.2f" %
                  (paymentDate,
                   self._fixedYearFracs[iFlow],
                   self._fixedFlows[iFlow],
                   self._fixedDfs[iFlow],
                   self._fixedFlowPVs[iFlow],
                   self._fixedTotalPV[iFlow]))

            iFlow += 1

##########################################################################

    def printFixedLegFlows(self):
        ''' Prints the fixed leg amounts without any valuation details. Shows
        the dates and sizes of the promised fixed leg flows. '''

        print("START DATE:", self._effectiveDate)
        print("MATURITY DATE:", self._maturityDate)
        print("COUPON (%):", self._fixedCoupon * 100)
        print("FIXED LEG FREQUENCY:", str(self._fixedFrequencyType))
        print("FIXED LEG DAY COUNT:", str(self._fixedDayCountType))

        if len(self._fixedFlows) == 0:
            print("Fixed Flows not calculated.")
            return

        header = "PAYMENT_DATE     YEAR_FRAC        FLOW"
        print(header)

        startIndex = 1

        iFlow = 0
        for paymentDate in self._adjustedFixedDates[startIndex:]:
            print("%15s %12.8f %12.2f" %
                  (paymentDate,
                   self._fixedYearFracs[iFlow],
                   self._fixedFlows[iFlow]))

            iFlow += 1

##########################################################################

    def printFloatLegPV(self):
        ''' Prints the floating leg dates, accrual factors, discount factors,
        forward libor rates, implied cash amounts, their present value and
        their cumulative PV using the last valuation performed. '''

        print("START DATE:", self._effectiveDate)
        print("MATURITY DATE:", self._maturityDate)
        print("SPREAD COUPON (%):", self._floatSpread * 100)
        print("FLOAT LEG FREQUENCY:", str(self._floatFrequencyType))
        print("FLOAT LEG DAY COUNT:", str(self._floatDayCountType))
        print("VALUATION DATE", self._valuationDate)

        if len(self._floatFlows) == 0:
            print("Floating Flows not calculated.")
            return

        if self._firstFixingRate is None:
            print("         *** FIRST FLOATING RATE PAYMENT IS IMPLIED ***")

        header = "PAYMENT_DATE     YEAR_FRAC    RATE(%)       FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        startIndex = self._floatStartIndex

        # By definition the discount factor is 1.0 on the valuation date

        print("%15s %10s %10s %12s %12.8f %12s %12s" %
              (self._valuationDate,
               "-",
               "-",
               "-",
               1.0,
               "-",
               "-"))

        iFlow = 0
        for paymentDate in self._adjustedFloatDates[startIndex:]:
            print("%15s %10.7f %10.5f %12.2f %12.8f %12.2f %12.2f" %
                  (paymentDate,
                   self._floatYearFracs[iFlow],
                   self._floatRates[iFlow]*100.0,
                   self._floatFlows[iFlow],
                   self._floatDfs[iFlow],
                   self._floatFlowPVs[iFlow],
                   self._floatTotalPV[iFlow]))

            iFlow += 1

##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START DATE", self._effectiveDate)
        s += labelToString("TERMINATION DATE", self._terminationDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("SWAP TYPE", self._swapType)
        s += labelToString("FIXED COUPON", self._fixedCoupon)
        s += labelToString("FLOAT SPREAD", self._floatSpread)
        s += labelToString("FIXED FREQUENCY", self._fixedFrequencyType)
        s += labelToString("FLOAT FREQUENCY", self._floatFrequencyType)
        s += labelToString("FIXED DAY COUNT", self._fixedDayCountType)
        s += labelToString("FLOAT DAY COUNT", self._floatDayCountType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUS DAY ADJUST", self._busDayAdjustType)
        s += labelToString("DATE GEN TYPE", self._dateGenRuleType)
        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################
