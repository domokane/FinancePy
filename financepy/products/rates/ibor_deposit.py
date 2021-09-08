##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.day_count import DayCount
from ...utils.day_count import DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types


###############################################################################


class IborDeposit:
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
                 start_date: Date,  # When the interest starts to accrue
                 maturity_date_or_tenor: (Date, str),  # Repayment of interest
                 deposit_rate: float,  # MM rate using simple interest
                 day_count_type: DayCountTypes,  # How year fraction is calculated
                 notional: float = 100.0,  # Amount borrowed
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,  # Holidays for maturity date
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.MODIFIED_FOLLOWING):
        """ Create a Libor deposit object which takes the start date when
        the amount of notional is borrowed, a maturity date or a tenor and the
        deposit rate. If a tenor is used then this is added to the start
        date and the calendar and business day adjustment method are applied if
        the maturity date fall on a holiday. Note that in order to calculate
        the start date you add the spot business days to the trade date
        which usually today. """

        check_argument_types(self.__init__, locals())

        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type

        if type(maturity_date_or_tenor) == Date:
            maturity_date = maturity_date_or_tenor
        else:
            maturity_date = start_date.add_tenor(maturity_date_or_tenor)

        calendar = Calendar(self._calendar_type)

        maturity_date = calendar.adjust(maturity_date,
                                        self._bus_day_adjust_type)
        if start_date > maturity_date:
            raise FinError("Start date cannot be after maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._deposit_rate = deposit_rate
        self._day_count_type = day_count_type
        self._notional = notional

    ###########################################################################

    def _maturity_df(self):
        """ Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date."""

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        df = 1.0 / (1.0 + acc_factor * self._deposit_rate)
        return df

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

        value = (1.0 + acc_factor * self._deposit_rate) * self._notional

        # Need to take into account spot days being zero so depo settling fwd
        value = value * df_maturity / df_settle

        return value

    ###########################################################################

    def print_flows(self,
                    valuation_date: Date):
        """ Print the date and size of the future repayment. """

        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        flow = (1.0 + acc_factor * self._deposit_rate) * self._notional
        print(self._maturity_date, flow)

    ###########################################################################

    def __repr__(self):
        """ Print the contractual details of the Libor deposit. """
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._start_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("DEPOSIT RATE", self._deposit_rate)
        s += label_to_string("DAY COUNT TYPE", self._day_count_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        s += label_to_string("BUS DAY ADJUST TYPE", self._bus_day_adjust_type)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
