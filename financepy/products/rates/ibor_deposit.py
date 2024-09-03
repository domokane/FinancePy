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
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types


###############################################################################


class IborDeposit:
    """An Ibor deposit is an agreement to borrow money interbank at the Ibor
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
    deposit contract."""

    def __init__(
        self,
        start_dt: Date,  # When the interest starts to accrue
        maturity_dt_or_tenor: (Date, str),  # Repayment of interest
        deposit_rate: float,  # MM rate using simple interest
        dc_type: DayCountTypes,  # How year fraction is calculated
        notional: float = 100.0,  # Amount borrowed
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,  # Maturity date
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.MODIFIED_FOLLOWING,
    ):
        """Create a Libor deposit object which takes the start date when
        the amount of notional is borrowed, a maturity date or a tenor and the
        deposit rate. If a tenor is used then this is added to the start
        date and the calendar and business day adjustment method are applied if
        the maturity date fall on a holiday. Note that in order to calculate
        the start date you add the spot business days to the trade date
        which usually today."""

        check_argument_types(self.__init__, locals())

        self.cal_type = cal_type
        self.bd_type = bd_type

        if type(maturity_dt_or_tenor) is Date:
            maturity_dt = maturity_dt_or_tenor
        else:
            maturity_dt = start_dt.add_tenor(maturity_dt_or_tenor)

        calendar = Calendar(self.cal_type)

        maturity_dt = calendar.adjust(maturity_dt, self.bd_type)

        if start_dt > maturity_dt:
            raise FinError("Start date cannot be after maturity date")

        self.start_dt = start_dt
        self.maturity_dt = maturity_dt
        self.deposit_rate = deposit_rate
        self.dc_type = dc_type
        self.notional = notional

    ###########################################################################

    def _maturity_df(self):
        """Returns the maturity date discount factor that would allow the
        Libor curve to reprice the contractual market deposit rate. Note that
        this is a forward discount factor that starts on settlement date."""

        dc = DayCount(self.dc_type)
        acc_factor = dc.year_frac(self.start_dt, self.maturity_dt)[0]
        df = 1.0 / (1.0 + acc_factor * self.deposit_rate)
        return df

    ###########################################################################

    def value(self, value_dt: Date, libor_curve):
        """Determine the value of an existing Libor Deposit contract given a
        valuation date and a Libor curve. This is simply the PV of the future
        repayment plus interest discounted on the current Libor curve."""

        if value_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        dc = DayCount(self.dc_type)
        acc_factor = dc.year_frac(self.start_dt, self.maturity_dt)[0]
        df_settle = libor_curve.df(self.start_dt)
        df_maturity = libor_curve.df(self.maturity_dt)

        value = (1.0 + acc_factor * self.deposit_rate) * self.notional

        # Need to take into account spot days being zero so depo settling fwd
        value = value * df_maturity / df_settle

        return value

    ###########################################################################

    def valuation_details(
        self,
        valuation_date: Date,
        discount_curve: DiscountCurve,
        index_curve: DiscountCurve = None,
    ):
        """
        A long-hand method that returns various details relevant to valuation
        in a dictionary. Slower than value(...) so should not be used when
        performance is important

        We want thre output dictionary to have  the same labels for different
        benchmarks (depos, fras, swaps) because we want to present them
        together so please do not stick new outputs into one of them only

        TODO: make a test of this
        """
        if valuation_date > self.maturity_dt:
            raise FinError("Start date after maturity date")

        dc = DayCount(self.dc_type)
        acc_factor = dc.year_frac(self.start_dt, self.maturity_dt)[0]
        df_settle = discount_curve.df(self.start_dt)
        df_maturity = discount_curve.df(self.maturity_dt)

        value = (1.0 + acc_factor * self.deposit_rate) * self.notional

        # Need to take into account spot days being zero so depo settling fwd
        value = (
            value * df_maturity / df_settle
        )  # VP: ??? this looks like a start_date - forward value? not spot value? why?

        out = {
            "type": type(self).__name__,
            "start_date": self.start_dt,
            "maturity_date": self.maturity_dt,
            "day_count_type": self.dc_type.name,
            "notional": self.notional,
            "contract_rate": self.deposit_rate,
            "market_rate": (df_settle / df_maturity - 1) / acc_factor,
            # for depo pvbp is actually negative: rates up, value down. but probably makes sense to report as positive, asif for a spot-starting fra
            "spot_pvbp": acc_factor * df_maturity,
            "fwd_pvbp": acc_factor * df_maturity / df_settle,
            "unit_value": value / self.notional,
            "value": value,
            # ignoring bus day adj type, calendar for now
        }
        return out

    ###########################################################################

    def print_flows(self, valuation_date: Date):
        """Print the date and size of the future repayment."""

        dc = DayCount(self.dc_type)
        acc_factor = dc.year_frac(self.start_dt, self.maturity_dt)[0]
        flow = (1.0 + acc_factor * self.deposit_rate) * self.notional
        print(self.maturity_dt, flow)

    ###########################################################################

    def __repr__(self):
        """Print the contractual details of the Libor deposit."""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.start_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("DEPOSIT RATE", self.deposit_rate)
        s += label_to_string("DAY COUNT TYPE", self.dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY ADJUST TYPE", self.bd_type)
        return s

    ###########################################################################

    def _print(self):
        print(self)


###############################################################################
