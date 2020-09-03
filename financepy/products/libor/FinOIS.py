##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinSchedule import FinSchedule
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinDate import FinDate, dailyWorkingDaySchedule
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinHelperFunctions import labelToString

###############################################################################
###############################################################################


class FinOIS(object):
    ''' Class for managing overnight index swaps. This is a swap contract in
    which a fixed payment leg is exchanged for a floating coupon leg. There
    is no exchange of par.

    The contract lasts from a start date to a specified maturity date.
    The fixed coupon is the OIS fixed rate which is set at contract initiation.

    The floating rate is not known until the end of each payment period. It is
    calculated at the end of the period as it is based on daily observations
    of the overnight index rate which are compounded according to a specific
    convention. Hence the OIS floating rate is determined by the history of the
    OIS rates.

    In its simplest form, there is just one fixed rate payment and one floating
    rate payment at contract maturity. However when the contract becomes longer
    than one year the floating and fixed payments become periodic.

    The value of the contract is the NPV of the two coupon streams.
    Discounting is done on a supplied OIS curve which is itself implied by
    the term structure of market OIS rates. '''

    def __init__(self,
                 startDate: FinDate,
                 maturityDate: FinDate,
                 fixedRate: float,
                 fixedFrequencyType: FinFrequencyTypes,
                 fixedDayCountType: FinDayCountTypes,
                 floatFrequencyType: FinFrequencyTypes = FinFrequencyTypes.ANNUAL,
                 floatDayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 payFixedLeg: bool = True,
                 notional: float = ONE_MILLION,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create OIS object. '''

        checkArgumentTypes(self.__init__, locals())

        if startDate > maturityDate:
            raise FinError("Start date after maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._payFixedLeg = payFixedLeg
        self._notional = notional

        self._fixedRate = fixedRate

        self._fixedFrequencyType = fixedFrequencyType
        self._floatFrequencyType = floatFrequencyType

        self._fixedDayCountType = fixedDayCountType
        self._floatDayCountType = floatDayCountType

        self._payFixedLeg = payFixedLeg

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

        # we only generate flows once we have a valuation date
        self._fixedFlows = []
        self._floatFlows = []

        # we only generate dates once we have a valuation date
        self._adjustedFixedDates = []
        self._adjustedFloatDates = []

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

        self.generatePaymentDates(valueDate)

        dayCounter = FinDayCount(self._fixedDayCountType)
        self._fixedFlows = []
        prevDt = valueDate
        for dt in self._adjustedFixedDates:
            flow = dayCounter.yearFrac(prevDt, dt)[0] * self._fixedRate
            self._fixedFlows.append(flow)
            prevDt = dt

    ###########################################################################

    def generateFloatLegFlows(self, valueDate, indexCurve):
        ''' Generate the payment amounts on floating leg implied by index
        curve '''
        self.generatePaymentDates(valueDate)

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

    def rate(self, oisDates, oisFixings):
        ''' Calculate the OIS rate implied rate from the history of fixings.
        '''

        if len(oisDates) != len(oisFixings):
            raise FinError("Dates and fixings must have same length.")

        prevDt = oisDates[0]
        cmpd = 1.0
        dayCounter = FinDayCount(self._dayCountType)

        for dt, fixing in zip(oisDates[1:], oisFixings[1:]):
            alpha = dayCounter.yearFrac(prevDt, dt)[0]
            cmpd *= (1.0 + fixing * alpha)

        alpha = dayCounter.yearFrac(oisDates[0], oisDates[-1])[0]
        rate = (cmpd - 1.0) / alpha
        return rate

    ###########################################################################

    def value(self, valueDate, discountCurve):
        ''' Value the interest rate swap on a value date given a single Libor
        discount curve. '''
        principal = 1.0

        fixedLegValue = self.fixedLegValue(valueDate,
                                           discountCurve,
                                           principal)

        floatLegValue = self.floatLegValue(valueDate,
                                           discountCurve,
                                           discountCurve,
                                           principal)

        value = fixedLegValue - floatLegValue

        if self._payFixedLeg is True:
            value = value * (-1.0)

        return value * self._notional

    ###########################################################################

    def fixedLegValue(self, valueDate, discountCurve, principal=0.0):

        self.generateFixedLegFlows(self._startDate)

        pv = 0.0
        df = 1.0

        for dt, flow in zip(self._adjustedFixedDates, self._fixedFlows):
            df = discountCurve.df(dt)
            pv += df * flow

        pv = pv + principal * df
        z0 = discountCurve.df(valueDate)
        pv = pv / z0
        return pv

    ###########################################################################

    def floatLegValue(self,
                      valueDate,
                      discountCurve,
                      indexCurve,
                      principal=0.0):
        ''' Value the floating leg with payments from an index curve and
        discounting based on a supplied discount curve. '''
        basis = FinDayCount(self._floatDayCountType)

        if self._floatFlows == []:
            self.generateFloatLegFlows(valueDate, indexCurve)

        dt1 = valueDate
        pv = 0.0

        for dt2, flow in zip(self._adjustedFloatDates, self._floatFlows):
            df = discountCurve.df(dt2)
            alpha = basis.yearFrac(dt1, dt2)[0]
            pv += df * flow * alpha
            dt1 = dt2

        pv = pv + df * principal
        return pv

    ###########################################################################

    def df(self,
           oisRate,
           startDate,
           endDate):
        ''' Calculate the OIS rate implied discount factor. '''

        df = 1.0
        compoundingDates = dailyWorkingDaySchedule(startDate, endDate)
        prevDt = startDate
        for dt in compoundingDates[1:]:
            dayCount = dt - prevDt
            df = df / (1.0 + oisRate * dayCount / 360.0)

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
