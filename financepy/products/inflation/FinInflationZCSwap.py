##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numpy import np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinInflationZCSwap(object):
    ''' Class for managing a zero coupon inflation swap. In this there are two
    legs each with a single payment. The fixed leg pays a predetermined amount
    and the inflation leg pays an amount linked to the realised CPI index. '''

    def __init__(self,
                 startDate: FinDate,  # The date the FRA starts to accrue
                 maturityDateOrTenor: (FinDate, str),  # End of the Libor rate period
                 fixedRate: float,  # The fixed leg growth rate
                 dayCountType: FinDayCountTypes,  # For interest period
                 notional: float = 100.0,
                 payFixedRate: bool = True,  # True if the FRA rate is being paid
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        ''' Create a Forward Rate Agreeement object. '''

        checkArgumentTypes(self.__init__, locals())

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = startDate.addTenor(maturityDateOrTenor)
            calendar = FinCalendar(self._calendarType)
            maturityDate = calendar.adjust(maturityDate,
                                           self._busDayAdjustType)

        if startDate > maturityDate:
            raise FinError("Settlement date after maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._fixedRate = fixedRate
        self._payFixedRate = payFixedRate
        self._dayCountType = dayCountType
        self._notional = notional

    ###########################################################################

    def value(self, valuationDate, liborCurve, realisedCPI ):
        ''' Determine mark to market value of a ZC inflation swap based on the
        actual CPI rate. '''

        df = liborCurve.df(self._maturityDate)
        years = (self._maturityDate - valuationDate) / gDaysInYear

        fixedLegGrowth = np.power(1.0 + self._inflationRate, years)
        inflationLegGrowth = np.power(1.0 + self._inflationRate, years)

        fixedLegValuePV = (fixedLegGrowth - 1.0) * df
        inflationLegValuePV = (inflationLegGrowth - 1.0) * df

        v = (fixedLegValuePV - inflationLegValuePV) * self._notional

        if self._payFixedRate:
            v *= -1.0

        return v

    ##########################################################################

    def maturityDeflator(self, liborCurve):

        dc = FinDayCount(self._dayCountType)
        df1 = liborCurve.df(self._startDate)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        df2 = df1 / (1.0 + accFactor * self._fraRate)
        return df2

    ###########################################################################

    def printFlows(self, valuationDate):
        ''' Determine the value of the Deposit given a Libor curve. '''

        flow_settle = self._notional
        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        flow_maturity = (1.0 + accFactor * self._fraRate) * self._notional

        if self._payFixedRate is True:
            print(self._startDate, -flow_settle)
            print(self._maturityDate, flow_maturity)
        else:
            print(self._startDate, flow_settle)
            print(self._maturityDate, -flow_maturity)

    ##########################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START ACCD DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("FIXED INFLATION RATE", self._fixedRate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("PAY FIXED RATE", self._payFixedRate)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        s += labelToString("BUS DAY ADJUST TYPE", self._busDayAdjustType)
        s += labelToString("CALENDAR", self._calendarType)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
