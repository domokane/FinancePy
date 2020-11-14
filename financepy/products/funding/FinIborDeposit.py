##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinCalendar import FinCalendar
from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinIborDeposit(object):
    ''' An Ibor deposit is an agreement to borrow money interbank at the Ibor
    fixing rate starting on the start date and repaid on the maturity date
    with the interest amount calculated according to a day count convention and
    dates calculated according to a calendar and business day adjustment rule.

    Care must be taken to calculate the correct start (settlement) date. Start
    with the trade (value) date which is typically today, we may need to add on
    a number of business days (spot days) to get to the settlement date. The
    maturity date is then calculated by adding on the deposit tenor/term to the
    settlement date and adjusting for weekends and holidays according to the
    calendar and adjustment type.

    Note that for over-night (ON) depos the settlement date is
    today with maturity in one business day. For tomorrow-next (TN) depos the
    settlement is in one business day with maturity on the following business
    day. For later maturity deposits, settlement is usually in 1-3 business
    days. The number of days depends on the currency and jurisdiction of the
    deposit contract. '''

    def __init__(self,
                 startDate: FinDate,  #  When the interest starts to accrue
                 maturityDateOrTenor: (FinDate, str),  # Repayment of interest
                 depositRate: float,   # MM rate using simple interest
                 dayCountType: FinDayCountTypes,  # How year fraction is calculated
                 notional: float = 100.0,   # Amount borrowed
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,  #  Holidays for maturity date
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        ''' Create a Libor deposit object which takes the start date when
        the amount of notional is borrowed, a maturity date or a tenor and the
        deposit rate. If a tenor is used then this is added to the start
        date and the calendar and business day adjustment method are applied if
        the maturity date fall on a holiday. Note that in order to calculate
        the start date you add the spot business days to the trade date
        which usually today. '''

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
            raise FinError("Start date cannot be after maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._depositRate = depositRate
        self._dayCountType = dayCountType
        self._notional = notional

    ###########################################################################

    def _maturityDf(self):
        ''' Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date.'''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        discountFactor = 1.0 / (1.0 + accFactor * self._depositRate)
        return discountFactor

    ###########################################################################

    def value(self,
              valuationDate: FinDate,
              liborCurve):
        ''' Determine the value of an existing Libor Deposit contract given a
        valuation date and a Libor curve. This is simply the PV of the future
        repayment plus interest discounted on the current Libor curve. '''

        if valuationDate > self._maturityDate:
            raise FinError("Start date after maturity date")

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        df_settle = liborCurve.df(self._startDate)
        df_maturity = liborCurve.df(self._maturityDate)

        value = (1.0 + accFactor * self._depositRate) * self._notional

        # Need to take into account spot days being zero so depo settling fwd
        value = value * df_maturity / df_settle

        return value

    ###########################################################################

    def printFlows(self,
                   valuationDate: FinDate):
        ''' Print the date and size of the future repayment. '''

        dc = FinDayCount(self._dayCountType)
        accFactor = dc.yearFrac(self._startDate, self._maturityDate)[0]
        flow = (1.0 + accFactor * self._depositRate) * self._notional
        print(self._maturityDate, flow)

    ###########################################################################

    def __repr__(self):
        ''' Print the contractual details of the Libor deposit. '''
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START DATE", self._startDate)
        s += labelToString("MATURITY DATE", self._maturityDate)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("DEPOSIT RATE", self._depositRate)
        s += labelToString("DAY COUNT TYPE", self._dayCountType)
        s += labelToString("CALENDAR", self._calendarType)
        s += labelToString("BUS DAY ADJUST TYPE", self._busDayAdjustType)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
