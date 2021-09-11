##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.error import FinError
from ...utils.date import Date
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve

###############################################################################


class IborFRA:
    """ Class for managing LIBOR forward rate agreements. A forward rate
    agreement is an agreement to exchange a fixed pre-agreed rate for a
    floating rate linked to LIBOR that is not known until some specified
    future fixing date. The FRA payment occurs on or soon after this date
    on the FRA settlement date. Typically the timing gap is two days.

    A FRA is used to hedge a Ibor quality loan or lend of some agreed
    notional amount. This period starts on the settlement date of the
    FRA and ends on the maturity date of the FRA. For example a 1x4 FRA
    relates to a Ibor starting in 1 month for a loan period ending in 4
    months. Hence it links to 3-month Ibor rate. The amount received by a
    payer of fixed rate at settlement is:

        acc(1,2) * (Ibor(1,2) - FRA RATE) / (1 + acc(0,1) x Ibor(0,1))

    So the value at time 0 is

        acc(1,2) * (FWD Ibor(1,2) - FRA RATE) x df(0,2)

    If the base date of the curve is before the value date then we
    forward adjust this amount to that value date. For simplicity I have
    assumed that the fixing date and the settlement date are the same date.
    This should be amended later.

    The valuation below incorporates a dual curve approach.
    """

    def __init__(self,
                 start_date: Date,  # The date the FRA starts to accrue
                 # End of the Ibor rate period
                 maturity_date_or_tenor: (Date, str),
                 fraRate: float,  # The fixed contractual FRA rate
                 day_count_type: DayCountTypes,  # For interest period
                 notional: float = 100.0,
                 payFixedRate: bool = True,  # True if the FRA rate is being paid
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.MODIFIED_FOLLOWING):
        """ Create a Forward Rate Agreement object. """

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
            raise FinError("Settlement date after maturity date")

        self._start_date = start_date
        self._maturity_date = maturity_date
        self._fraRate = fraRate
        self._payFixedRate = payFixedRate
        self._day_count_type = day_count_type
        self._notional = notional

    ###########################################################################

    def value(self,
              valuation_date: Date,
              discount_curve: DiscountCurve,
              index_curve: DiscountCurve = None):
        """ Determine mark to market value of a FRA contract based on the
        market FRA rate. We allow the pricing to have a different curve for
        the Libor index and the discounting of promised cash flows. """

        if index_curve is None:
            index_curve = discount_curve

        # Get the Libor index from the index curve
        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        dfIndex1 = index_curve.df(self._start_date)
        dfIndex2 = index_curve.df(self._maturity_date)
        liborFwd = (dfIndex1 / dfIndex2 - 1.0) / acc_factor

        # Get the discount factor from a discount curve
        dfDiscount2 = discount_curve.df(self._maturity_date)

        v = acc_factor * (liborFwd - self._fraRate) * dfDiscount2

        # Forward value the FRA to the value date
        df_to_valuation_date = discount_curve.df(valuation_date)
        v = v * self._notional / df_to_valuation_date

        if self._payFixedRate is True:
            v *= -1.0
        return v

    ##########################################################################

    def maturity_df(self, index_curve):
        """ Determine the maturity date index discount factor needed to refit
        the market FRA rate. In a dual-curve world, this is not the discount
        rate discount factor but the index curve discount factor. """

        dc = DayCount(self._day_count_type)
        df1 = index_curve.df(self._start_date)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        df2 = df1 / (1.0 + acc_factor * self._fraRate)
        return df2

    ###########################################################################

    def print_flows(self, valuation_date):
        """ Determine the value of the Deposit given a Ibor curve. """

        flow_settle = self._notional
        dc = DayCount(self._day_count_type)
        acc_factor = dc.year_frac(self._start_date, self._maturity_date)[0]
        flow_maturity = (1.0 + acc_factor * self._fraRate) * self._notional

        if self._payFixedRate is True:
            print(self._start_date, -flow_settle)
            print(self._maturity_date, flow_maturity)
        else:
            print(self._start_date, flow_settle)
            print(self._maturity_date, -flow_maturity)

    ##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START ACCD DATE", self._start_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("FRA RATE", self._fraRate)
        s += label_to_string("NOTIONAL", self._notional)
        s += label_to_string("PAY FIXED RATE", self._payFixedRate)
        s += label_to_string("DAY COUNT TYPE", self._day_count_type)
        s += label_to_string("BUS DAY ADJUST TYPE", self._bus_day_adjust_type)
        s += label_to_string("CALENDAR", self._calendar_type)
        return s

    ###########################################################################

    def _print(self):
        print(self)

###############################################################################
