##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinCalendar import FinCalendarTypes,  FinDateGenRuleTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinHelperFunctions import labelToString

##########################################################################


class FinLiborSwap(object):
    ''' Class for managing an interest rate swap contract. '''

    def __init__(self,
                 startDate,  # This is typically T+2 on a new swap
                 maturityDateOrTenor,
                 fixedCoupon,
                 fixedFreqType,
                 fixedDayCountType,
                 notional=100.0,
                 floatSpread=0.0,
                 floatFreqType=FinFrequencyTypes.QUARTERLY,
                 floatDayCountType=FinDayCountTypes.THIRTY_360,
                 payFixedFlag=True,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):
        ''' Create an interest rate swap contract. '''

        if type(startDate) != FinDate:
            raise ValueError("Settlement date must be a FinDate.")

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = startDate.addTenor(maturityDateOrTenor)

        if startDate > maturityDate:
            raise ValueError("Start date after maturity date")

        if fixedDayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Fixed Day Count Rule type " +
                str(fixedDayCountType))

        if floatDayCountType not in FinDayCountTypes:
            raise ValueError(
                "Unknown Float Day Count Rule type " +
                str(floatDayCountType))

        if fixedFreqType not in FinFrequencyTypes:
            raise ValueError(
                "Unknown Fixed Frequency type " +
                str(fixedFreqType))

        if floatFreqType not in FinFrequencyTypes:
            raise ValueError(
                "Unknown Float Frequency type " +
                str(fixedFreqType))

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinBusDayAdjustTypes:
            raise ValueError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise ValueError(
                "Unknown Date Gen Rule type " +
                str(dateGenRuleType))

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._notional = notional
        self._payFixedLeg = payFixedFlag

        self._fixedCoupon = fixedCoupon
        self._floatSpread = floatSpread

        self._fixedFrequencyType = fixedFreqType
        self._floatFrequencyType = floatFreqType

        self._fixedDayCountType = fixedDayCountType
        self._floatDayCountType = floatDayCountType

        self._payFixedFlag = payFixedFlag  ## IS THIS A DUPLICATE OF ABOVE

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        # These are generated immediately as they are for the entire
        # life of the swap. Given a valuation date we can determine
        # which cash flows are in the future and value the swap
        self._generateFixedLegPaymentDates()
        self._generateFloatLegPaymentDates()

        # Need to know latest payment date for bootstrap
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

        if self._payFixedLeg is True:
            value = value * (-1.0)

        return value

##########################################################################

    def _generateFixedLegPaymentDates(self):
        ''' Generate the fixed leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFixedDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._fixedFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType).generate()

##########################################################################

    def _generateFloatLegPaymentDates(self):
        ''' Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date'''
        self._adjustedFloatDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._floatFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType).generate()

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

    def parCoupon(self, valuationDate, discountCurve):
        ''' Calculate the fixed leg coupon that makes the swap worth zero.
        If the valuation date is before the swap payments start then this
        is the forward swap rate as it starts in the future. The swap rate
        is then a forward swap rate and so we use a forward discount
        factor. If the swap fixed leg has begun then we have a spot
        starting swap. '''

        pv01 = self.pv01(valuationDate, discountCurve)

        if valuationDate < self._startDate:
            df0 = discountCurve.df(self._startDate)
        else:
            df0 = discountCurve.df(valuationDate)

        dfT = discountCurve.df(self._maturityDate)
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

        basis = FinDayCount(self._fixedDayCountType)

        ''' The swap may have started in the past but we can only value
        payments that have occurred after the valuation date. '''
        startIndex = 0
        while self._adjustedFixedDates[startIndex] < valuationDate:
            startIndex += 1

        ''' If the swap has yet to settle then we do not include the
        start date of the swap as a coupon payment date. '''
        if valuationDate <= self._startDate:
            startIndex = 1

        self._fixedStartIndex = startIndex

        ''' Now PV fixed leg flows. '''
        pv = 0.0
        prevDt = self._adjustedFixedDates[startIndex - 1]

        for nextDt in self._adjustedFixedDates[startIndex:]:
            alpha = basis.yearFrac(prevDt, nextDt)
            df_discount = discountCurve.df(nextDt)
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

        z0 = discountCurve.df(valuationDate)
        pv = pv / z0
        return pv

##########################################################################

    def floatLegValue(self,
                      valuationDate,  # IS THIS THE SETTLEMENT DATE ????
                      discountCurve,
                      indexCurve,
                      firstFixingRate=None,
                      principal=0.0):
        ''' Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. '''

        self._valuationDate = valuationDate
        self._floatYearFracs = []
        self._floatFlows = []
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

        ''' The first floating payment is usually already fixed so is
        not implied by the index curve. '''
        prevDt = self._adjustedFloatDates[startIndex - 1]
        nextDt = self._adjustedFloatDates[startIndex]
        alpha = basis.yearFrac(prevDt, nextDt)
        df1_index = 1.0  # Cannot be previous date as that has past
        df2_index = indexCurve.df(nextDt)

        if self._firstFixingRate is None:
            libor = (df1_index / df2_index - 1.0) / alpha
            flow = libor * alpha * self._notional
        else:
            flow = self._firstFixingRate * alpha * self._notional

        df_discount = discountCurve.df(nextDt)
        pv = flow * df_discount

        self._floatYearFracs.append(alpha)
        self._floatFlows.append(flow)
        self._floatDfs.append(df_discount)
        self._floatFlowPVs.append(flow * df_discount)
        self._floatTotalPV.append(pv)

        prevDt = nextDt
        df1_index = indexCurve.df(prevDt)

        for nextDt in self._adjustedFloatDates[startIndex + 1:]:
            alpha = basis.yearFrac(prevDt, nextDt)
            df2_index = indexCurve.df(nextDt)
            # The accrual factors cancel
            flow = (df1_index / df2_index - 1.0) * self._notional
            df_discount = discountCurve.df(nextDt)
            pv += flow * df_discount
            df1_index = df2_index
            prevDt = nextDt

            self._floatFlows.append(flow)
            self._floatYearFracs.append(alpha)
            self._floatDfs.append(df_discount)
            self._floatFlowPVs.append(flow * df_discount)
            self._floatTotalPV.append(pv)

        flow = principal * self._notional
        pv = pv + flow * df_discount
        self._floatFlows[-1] += flow
        self._floatFlowPVs[-1] += flow * df_discount
        self._floatTotalPV[-1] = pv

        z0 = discountCurve.df(valuationDate)
        pv = pv / z0

        return pv

##########################################################################

    def printFixedLeg(self):
        ''' Prints the fixed leg amounts. '''

        print("START DATE:", self._startDate)
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

        startIndex = self._fixedStartIndex

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

    def printFloatLeg(self):
        ''' Prints the floating leg amounts. '''

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

        header = "PAYMENT_DATE     YEAR_FRAC        FLOW         DF"
        header += "         DF*FLOW       CUM_PV"
        print(header)

        startIndex = self._floatStartIndex

        iFlow = 0
        for paymentDate in self._adjustedFloatDates[startIndex:]:
            print("%15s %10.7f %12.2f %12.8f %12.2f %12.2f" %
                  (paymentDate,
                   self._floatYearFracs[iFlow],
                   self._floatFlows[iFlow],
                   self._floatDfs[iFlow],
                   self._floatFlowPVs[iFlow],
                   self._floatTotalPV[iFlow]))

            iFlow += 1

##########################################################################

    def __repr__(self):
        s = labelToString("START DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("PAY FIXED LEG", self._payFixedLeg)
        s += labelToString("FIXED COUPON", self._fixedCoupon)
        s += labelToString("FLOAT SPREAD", self._floatSpread)
        s += labelToString("FIXED FREQUENCY", self._fixedFrequencyType)
        s += labelToString("FLOAT FREQUENCY", self._floatFrequencyType)
        s += labelToString("FIXED DAY COUNT", self._fixedDayCountType)
        s += labelToString("FLOAT DAY COUNT", self._floatDayCountType)
        s += labelToString("PAY FIXED FLAG", self._payFixedFlag)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUS DAY ADJUST", self._busDayAdjustType)
        s += labelToString("DATE GEN TYPE", self._dateGenRuleType)
        return s

###############################################################################

    def print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################
