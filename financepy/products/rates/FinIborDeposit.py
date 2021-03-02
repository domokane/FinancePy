##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.Date import Date
from ...utils.FinError import FinError
from ...utils.Calendar import Calendar
from ...utils.Calendar import FinCalendarTypes
from ...utils.Calendar import FinBusDayAdjustTypes
from ...utils.DayCount import DayCount
from ...utils.DayCount import FinDayCountTypes
from ...utils.FinHelperFunctions import labelToString, checkArgumentTypes

###############################################################################


class FinIborDeposit(object):
    """ An Ibor deposit is an agreement to borrow money interbank at the Ibor
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
    deposit contract. """

    def __init__(self,
                 start_date: Date,  #  When the interest starts to accrue
                 maturity_date_or_tenor: (Date, str),  # Repayment of interest
                 depositRate: float,  # MM rate using simple interest
                 day_count_type: FinDayCountTypes,  # How year fraction is calculated
                 notional: float = 100.0,  # Amount borrowed
                 calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND,  #  Holidays for maturity date
                 bus_day_adjust_type: FinBusDayAdjustTypes = FinBusDayAdjustTypes.MODIFIED_FOLLOWING):
        """ Create a Libor deposit object which takes the start date when
        the amount of notional is borrowed, a maturity date or a tenor and the
        deposit rate. If a tenor is used then this is added to the start
        date and the calendar and business day adjustment method are applied if
        the maturity date fall on a holiday. Note that in order to calculate
        the start date you add the spot business days to the trade date
        which usually today. """

        checkArgumentTypes(self.__init__, locals())

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            maturity_date = start_date.addTenor(maturity_date_or_tenor)

        calendar = Calendar(self._calendar_type)

        maturity_date = calendar.adjust(maturity_date,
                                       self._bus_day_adjust_type)
        if start_date > maturity_date:
            raise FinError("Start date cannot be after maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._depositRate = depositRate
        self._day_count_type = day_count_type
        self._notional = notional

    ###########################################################################

    def _maturityDf(self):
        """ Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date."""

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        discountFactor = 1.0 / (1.0 + acc_factor * self._depositRate)
        return discountFactor

    ###########################################################################

    def value(self,
              valuation_date: Date,
              libor_curve):
        """ Determine the value of an existing Libor Deposit contract given a
        valuation date and a Libor curve. This is simply the PV of the future
        repayment plus interest discounted on the current Libor curve. """

        if valuation_date > self._maturity_date:
            raise FinError("Start date after maturity date")

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        df_settle = libor_curve.df(self._start_date)
        df_maturity = libor_curve.df(self._maturity_date)

        value = (1.0 + acc_factor * self._depositRate) * self._notional

        # Need to take into account spot days being zero so depo settling fwd
        value = value * df_maturity / df_settle

        return value

    ###########################################################################

    def printFlows(self,
                   valuation_date: Date):
        """ Print the date and size of the future repayment. """

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        flow = (1.0 + acc_factor * self._depositRate) * self._notional
        print(self._maturity_date, flow)

    ###########################################################################

    def __repr__(self):
        """ Print the contractual details of the Libor deposit. """
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("START DATE", self._start_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("NOTIONAL", self._notional)
        s += labelToString("DEPOSIT RATE", self._depositRate)
        s += labelToString("DAY COUNT TYPE", self._day_count_type)
        s += labelToString("CALENDAR", self._calendar_type)
        s += labelToString("BUS DAY ADJUST TYPE", self._bus_day_adjust_type)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
