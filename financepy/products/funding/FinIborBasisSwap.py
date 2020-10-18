##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes, FinDateGenRuleTypes
from ...finutils.FinCalendar import FinCalendar, FinBusDayAdjustTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalTypes import FinSwapTypes

##########################################################################


class FinIborBasisSwap(object):
    ''' Class for managing an interest rate basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a 
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.
    
    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on a supplied discount curve which is separate from the curves from
    which the implied index rates are extracted. '''
    
    def __init__(self,
                 startDate: FinDate,  # Date interest starts to accrue
                 terminationDateOrTenor: (FinDate, str),  # Date contract ends
                 payFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 payDayCountType: FinDayCountTypes  = FinDayCountTypes.THIRTY_E_360,
                 recFreqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 recDayCountType: FinDayCountTypes = FinDayCountTypes.THIRTY_E_360,
                 basisSwapSpread: float = 0.0,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating 
        legs and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated. '''

        checkArgumentTypes(self.__init__, locals())

        if type(terminationDateOrTenor) == FinDate:
            self._terminationDate = terminationDateOrTenor
        else:
            self._terminationDate = startDate.addTenor(terminationDateOrTenor)

        calendar = FinCalendar(calendarType)
        self._maturityDate = calendar.adjust(self._terminationDate,
                                             busDayAdjustType)

        if startDate > self._maturityDate:
            raise FinError("Start date after maturity date")

        self._startDate = startDate
        self._notional = notional
        self._basisSwapSpread = basisSwapSpread
        self._payFreqType = payFreqType
        self._recFreqType = recFreqType
        self._payDayCountType = payDayCountType
        self._recDayCountType = recDayCountType
        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        self._payFloatDates = self._generateFloatLegDates(self._payFreqType)
        self._recFloatDates = self._generateFloatLegDates(self._recFreqType)

        self._adjustedMaturityDate = self._adjustedFixedDates[-1]

        self._payFloatYearFracs = []
        self._recFloatYearFracs = []

        self._payFloatFlows = []
        self._recFloatFlows = []

        self._payFloatFlowPVs = []
        self._recFloatFlowPVs = []

        self._payFirstFixingRate = None
        self._recFirstFixingRate = None

        self._valuationDate = None

##########################################################################

    def _generatePayFloatLegPaymentDates(self, frequencyType):
        ''' Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        
        floatDates = FinSchedule(self._startDate,
                                 self._terminationDate,
                                 frequencyType,
                                 self._calendarType,
                                 self._busDayAdjustType,
                                 self._dateGenRuleType)._generate()

        return floatDates

##########################################################################

    def value(self,
              valuationDate,
              discountCurve,
              payIndexCurve,
              recIndexCurve,
              payFirstFixingRate=None,
              recFirstFixingRate=None,
              principal=0.0):

        ''' Value the LIBOR basis swap on a value date given a single Ibor
        discount curve and each of the index curves for the two floating legs
        of the swap. '''

        payFloatLegValue = self.floatLegValue(valuationDate,
                                              discountCurve,
                                              payIndexCurve,
                                              payFirstFixingRate,
                                              principal)

        recFloatLegValue = self.floatLegValue(valuationDate,
                                              discountCurve,
                                              recIndexCurve,
                                              recFirstFixingRate,
                                              principal)

        value = recFloatLegValue - payFloatLegValue
        return value

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
        if valuationDate <= self._startDate:
            startIndex = 1

        self._floatStartIndex = startIndex

        # Forward price to settlement date (if valuation is settlement date)
        self._dfValuationDate = discountCurve.df(valuationDate)

        ''' The first floating payment is usually already fixed so is
        not implied by the index curve. '''
        prevDt = self._adjustedFloatDates[startIndex - 1]
        nextDt = self._adjustedFloatDates[startIndex]
        alpha = basis.yearFrac(prevDt, nextDt)[0]
        df1_index = indexCurve.df(self._startDate)  # Cannot be pcd as has past
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

    def printFloatLegPV(self):
        ''' Prints the floating leg dates, accrual factors, discount factors,
        forward libor rates, implied cash amounts, their present value and
        their cumulative PV using the last valuation performed. '''

        print("START DATE:", self._startDate)
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
        s += labelToString("START DATE", self._startDate)
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
