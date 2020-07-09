##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinLiborDeposit(object):
    ''' A Libor deposit is an agreement to borrow money interbank at the Libor
    fixing rate starting on the settlement date and repaid on the maturity date
    with the interest amount calculated according to a day count convention and
    dates calculated according to a calendar and business day adjustment rule.
    '''

    def __init__(self,
                 settlementDate: FinDate,
                 maturityDateOrTenor: Union[FinDate, str],
                 depositRate: float,
                 dayCountType: FinDayCountTypes,
                 notional: float = 100.0,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        ''' Create a Libor deposit object which takes the settlement date when
        the amount of notional is borrowed, the deposit rate, the day count
        convention used to calculate the interest paid and a calendar and a
        business day adjustment method if dates fall on holidays. '''

        checkArgumentTypes(self.__init__, locals())

        self._calendarType = calendarType
        self._busDayAdjustType = busDayAdjustType

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = settlementDate.addTenor(maturityDateOrTenor)
            calendar = FinCalendar(self._calendarType)
            maturityDate = calendar.adjust(maturityDate,
                                           self._busDayAdjustType)

        if settlementDate > maturityDate:
            raise FinError("Settlement date after maturity date")

        self._settlementDate = settlementDate
        self._maturityDate = maturityDate
        self._depositRate = depositRate
        self._dayCountType = dayCountType
        self._notional = notional

    ###########################################################################

    def maturityDf(self):
        ''' Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date.'''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        discountFactor = 1.0 / (1.0 + accFactor * self._depositRate)
        return discountFactor

    ###########################################################################

    def value(self, valuationDate, liborCurve):
        ''' Determine the value of an existing Libor Deposit contract given a
        valuation date and a Libor curve. This is simply the PV of the future
        repayment plus interest discounted on the current Libor curve. '''

        if valuationDate > self._maturityDate:
            raise FinError("Start date after maturity date")

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        df = liborCurve.df(self._maturityDate)
        value = (1.0 + accFactor * self._depositRate) * df * self._notional

        df_valueDate = liborCurve.df(valuationDate)
        value = value / df_valueDate

        return value

    ###########################################################################

    def printFlows(self, valuationDate):
        ''' Print the date and size of the future repayment. '''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._settlementDate, self._maturityDate)
        flow = (1.0 + accFactor * self._depositRate) * self._notional
        print(self._maturityDate, flow)

    ###########################################################################

    def __repr__(self):
        ''' Print the contractual details of the Libor deposit. '''
        s = labelToString("SETTLEMENT DATE", self._settlementDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("DEPOSIT RATE", self._depositRate)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUS DAY ADJUST TYPE", self._busDayAdjustType)
        return s

    ###########################################################################

    def print(self):
        print(self)

###############################################################################
