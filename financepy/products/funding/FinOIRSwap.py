##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinDate import FinDate
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinGlobalTypes import FinSwapTypes

###############################################################################

from enum import Enum

class FinCompoundingTypes(Enum):
     COMPOUNDED = 1
     OVERNIGHT_COMPOUNDED_ANNUAL_RATE = 2
     AVERAGED = 3
     AVERAGED_DAILY = 4

###############################################################################


class FinOIRSwap(object):
    ''' Class for managing overnight index rate swaps (OIS) and Fed Fund swaps. 
    This is a contract in which a fixed payment leg is exchanged for a payment
    which pays the rolled-up overnight index rate. There is no exchange of par.
    The contract is entered into at zero initial cost.

    The contract lasts from a start date to a specified maturity date.
    The fixed coupon is the OIS fixed rate for the corresponding tenor which is
    set at contract initiation.

    The floating rate is not known fully until the end of each payment period.
    It's calculated at the contract maturity and is based on daily observations
    of the overnight index rate which are compounded according to a specific
    convention. Hence the OIS floating rate is determined by the history of the
    OIS rates.

    In its simplest form, there is just one fixed rate payment and one floating
    rate payment at contract maturity. However when the contract becomes longer
    than one year the floating and fixed payments become periodic, usually with
    annual exchanges of cash.

    The value of the contract is the NPV of the two coupon streams. Discounting
    is done on the OIS curve which is itself implied by the term structure of
    market OIS rates. '''

    def __init__(self,
                 startDate: FinDate,
                 maturityDate: FinDate,
                 swapType: FinSwapTypes,
                 fixedRate: float,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 floatFrequencyType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 notional: float = ONE_MILLION,
                 compoundMethod: FinCompoundingTypes = FinCompoundingTypes.COMPOUNDED,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create OIS object. '''

        checkArgumentTypes(self.__init__, locals())

        if startDate > maturityDate:
            raise FinError("Start date after maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._swapType = swapType
        self._notional = notional

        self._fixedRate = fixedRate

        self._fixedFrequencyType = fixedFrequencyType
        self._floatFrequencyType = floatFrequencyType

        self._fixedDayCountType = fixedDayCountType
        self._floatDayCountType = floatDayCountType

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        # we only generate flows once we have a valuation date
        self._fixedFlows = []
        self._floatFlows = []

        # There is only one payment date at the contract maturity
        self._adjustedFixedDates = [self._maturityDate]
        self._adjustedFloatDates = [self._maturityDate]

    ###########################################################################

    def generatePaymentDates(self, valueDate):

        self._adjustedFixedDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._fixedFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType)._generate()

        self._adjustedFloatDates = FinSchedule(
            self._startDate,
            self._maturityDate,
            self._floatFrequencyType,
            self._calendarType,
            self._busDayAdjustType,
            self._dateGenRuleType)._generate()

    ###########################################################################

    def generateFixedLegFlows(self, valueDate):

        dayCounter = FinDayCount(self._fixedDayCountType)
        yearFrac = dayCounter.yearFrac(self._startDate, self._maturityDate)[0]
        flow = yearFrac * self._fixedRate
        self._fixedFlows.append(flow)

    ###########################################################################

    def generateFloatLegFlows(self, valueDate, indexCurve):
        ''' Generate the payment amounts on floating leg implied by index
        curve. This requires us to calculate the daily rate and roll it up. '''

        dayCounter = FinDayCount(self._floatDayCountType)

        dt1 = valueDate
        df1 = indexCurve.df(dt1)

        self._floatFlows = []
        prevDt = valueDate

        for dt in self._adjustedFloatDates[1:]:
            alpha = dayCounter.yearFrac(prevDt, dt)[0]
            df2 = indexCurve.df(dt)
            flow = (df1 / df2 - 1.0) / alpha
            self._floatFlows.append(flow)
            prevDt = dt

    ###########################################################################

    def fixedLegValue(self, valueDate, discountCurve):
        ''' Calculate the OIS growth implied rate from the fixed rate which is
        compounded daily from the start of the overnight swap until maturity
        and then discounted back to the value date. If there are multiple pay-
        ments then we grow the '''

        dayCounter = FinDayCount(self._fixedDayCountType)
        dt = self._startDate
        value = 0.0

        for cpnDate in self._adjustedFixedDates:

            cmpd = 1.0
    
            while dt <= cpnDate:
                prevDt = dt
                dt = dt.addWeekDays(1)
                alpha = dayCounter.yearFrac(prevDt, dt)[0]
                cmpd *= (1.0 + self._fixedRate * alpha)

            df = discountCurve.df(cpnDate)

            value += self._notional * cmpd * df

        return value

    ###########################################################################

    def floatLegValue(self, valueDate, prevDates, prevOvernightRates, 
                      oisCurve, discountCurve):
        ''' Calculate the value of the OIS floating leg at a value date in the
        trade when some of the past OIS rates will have set. '''

        numPrevDates = len(prevDates)
        if numPrevDates != len(prevOvernightRates):
            raise FinError("Number previous dates not equal to number of rates.")

        dayCounter = FinDayCount(self._floatDayCountType)

        cmpd = 1.0

        prevDt = self._startDate
        nextDt = prevDt.addWeekDays(1)

        if numPrevDates > 0:
 
            if prevDates[0] != nextDt:
                raise FinError("Prev dates do not start just after OIS start.")
    
            if prevDates[-1] != valueDate.addWeekDays(-1):
                raise FinError("Prev dates do not end just before value date.")

            for dt in prevDates[1:]:
                dt = dt.addWeekDays(1)
                alpha = dayCounter.yearFrac(prevDt, dt)[0]
                cmpd *= (1.0 + self._fixedRate * alpha)
                prevDt = dt

        dt = valueDate.addWeekDays(1)

        # Now consider growth from value date
        while dt <= self._maturityDate:
            dt = dt.addWeekDays(1)
            alpha = dayCounter.yearFrac(prevDt, dt)[0]
            cmpd *= (1.0 + self._fixedRate * alpha)
            prevDt = dt

        df = discountCurve.df(self._maturityDate)
        value = self._notional * cmpd * df
        return value

    ###########################################################################

    def value(self, valueDate,
              prevDates, prevRates, 
              discountCurve, oisCurve):
        ''' Value the OIS floating leg on a value date given the history of
        previous fixings dates and fixing rates and the current term structure
        of discount rates and ois rates. '''

        fixedLegValue = self.fixedLegValue(valueDate, discountCurve)

        floatLegValue = self.floatLegValue(valueDate, prevDates, prevRates,
                                           oisCurve, discountCurve)

        value = fixedLegValue - floatLegValue

        if self._swapType is FinSwapTypes.PAYER:
            value = value * (-1.0)

        return value * self._notional

    ###########################################################################

    def df(self, oisRate: float):
        ''' Calculate the OIS rate implied discount factor. '''

        df = 1.0
        prevDt = self._startDate
        dt = prevDt

        while dt <= self._maturityDate:
            dt = dt.addWeekDays(1)
            dayCount = dt - prevDt
            df = df / (1.0 + oisRate * dayCount / 360.0)
            prevDt = dt            

        return df

    ###########################################################################

    def printFlows(self, valueDate, indexCurve):
        ''' Print the dates and cash flows on the OIS. '''
        print("StartDate:", self._startDate)
        print("MaturityDate:", self._maturityDate)
        print("OISFixedRate:", self._fixedRate)
        print("PayFixedLeg:", self._payFixedLeg)
        print("Notional:", self._notional)
        print("FixedDayCountType:", self._fixedDayCountType)
        print("FixedFrequencyType:", self._fixedFrequencyType)
        print("FloatDayCountType:", self._floatDayCountType)
        print("FloatFrequencyType:", self._floatFrequencyType)

        self.generateFixedLegFlows(valueDate)
        self.generateFloatLegFlows(valueDate, indexCurve)

        print("Fixed Leg Flows")
        for dt in self._adjustedFixedDates:
            print(dt)

        print("Floating Leg Flows")
        for dt in self._adjustedFloatDates:
            print(dt)

    ###########################################################################

    def __repr__(self):
        ''' Print the dates and cash flows on the OIS. '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("StartDate:", self._startDate)
        s += labelToString("MaturityDate:", self._maturityDate)
        s += labelToString("OISFixedRate:", self._fixedRate)
        s += labelToString("PayFixedLeg:", self._payFixedLeg)
        s += labelToString("Notional:", self._notional)
        s += labelToString("FixedDayCountType:", self._fixedDayCountType)
        s += labelToString("FixedFrequencyType:", self._fixedFrequencyType)
        s += labelToString("FloatDayCountType:", self._floatDayCountType)
        s += labelToString("FloatFrequencyType:", self._floatFrequencyType)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
