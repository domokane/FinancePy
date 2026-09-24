##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...utils.global_types import OptionTypes

from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...utils.schedule import Schedule
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve
from ...utils.check_values import check_curve_dt
from ...utils.helpers import option_years

from ...models.black_scholes_analytic import european_value
from ...models.black_scholes import BlackScholes
from ...models.model import Model

########################################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
########################################################################################


class EquityForwardStartOption(EquityOption):
    """An EquityForwardStartOption is an option whose strike is set equal to the
    stock price on the start date and which expires at a later date. When a
    frequency is given, the period between the start date and the final expiry
    date is divided into a sequence of consecutive forward-start options, each
    one starting when the previous one expires with its strike set to the stock
    price at that time. With a single period this is the classic forward-start
    option of Rubinstein (1991); with several periods it is a cliquet of
    at-the-money forward-start options whose first strike is set on the start
    date rather than on the valuation date."""

    def __init__(
        self,
        start_dt: Date,
        final_expiry_dt: Date,
        opt_type: OptionTypes,
        freq_type: FrequencyTypes,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create the EquityForwardStartOption by passing in the start date on
        which the first strike is set, the final expiry date and whether it is a
        call or a put. The frequency determines the reset dates between the two
        dates; the calendar and business day conventions adjust those dates."""

        check_argument_types(self.__init__, locals())

        if opt_type != OptionTypes.EUROPEAN_CALL and opt_type != OptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown OPTION_TYPE" + str(opt_type))

        if final_expiry_dt <= start_dt:
            raise FinError("Expiry date must be after start date")

        self.start_dt = start_dt
        self.final_expiry_dt = final_expiry_dt
        self.opt_type = opt_type
        self.freq_type = freq_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        self.v_options = None
        self.dfs = None
        self.actual_dts = None

        self.expiry_dts = Schedule(
            self.start_dt,
            self.final_expiry_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        ).adjusted_dts

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Value the forward-start option as a sequence of at-the-money
        forward-start options using the Black-Scholes model. The option that
        starts on date t1 and expires on date t2 is worth
        S * Q(0, t1) * BS(1, t2 - t1, 1, r(t1, t2), q(t1, t2), vol) where
        Q(0, t1) is the dividend deflator to the strike setting date and r and
        q are the forward rates over the life of the option."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if value_dt > self.final_expiry_dt:
            raise FinError("Value date after final expiry date.")

        if value_dt > self.start_dt:
            raise FinError(
                "Value date after start date: the strike has been set and the "
                "option should be valued as a vanilla option."
            )

        if not isinstance(model, BlackScholes):
            raise FinError("Unknown Model Type")

        s0 = stock_price
        fwd_vol = max(model.volatility, 1e-6)

        v_total = 0.0
        self.v_options = []
        self.dfs = []
        self.actual_dts = []

        dt_prev = self.start_dt

        for dt in self.expiry_dts:

            if dt <= dt_prev:
                continue

            t_vol = option_years(dt_prev, dt)
            df_end = discount_curve.df(dt)

            # The deflator is out to the strike setting date
            dq_start = dividend_curve.df(dt_prev)

            # The rates are the forward rates over the life of the option
            fwd_r = discount_curve.fwd_zero_rate_cc(dt_prev, dt)
            fwd_q = dividend_curve.fwd_zero_rate_cc(dt_prev, dt)

            v = european_value(1.0, t_vol, 1.0, fwd_r, fwd_q, fwd_vol, self.opt_type.value)

            v_fwd_opt = s0 * dq_start * v
            v_total += v_fwd_opt

            self.dfs.append(df_end)
            self.v_options.append(v_fwd_opt)
            self.actual_dts.append(dt)

            dt_prev = dt

        return v_total

    ###########################################################################

    def print_payments(self):
        if self.v_options is None:
            raise FinError("Options not created yet.")
        num_options = len(self.v_options)
        for i in range(0, num_options):
            print(self.actual_dts[i], self.dfs[i], self.v_options[i])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START_DATE", self.start_dt)
        s += label_to_string("FINAL EXPIRY DATE", self.final_expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("FREQUENCY TYPE", self.freq_type)
        s += label_to_string("CALENDAR TYPE", self.cal_type)
        s += label_to_string("BUS_DAY_ADJUST", self.bd_type)
        s += label_to_string("DATE GEN RULE TYPE", self.dg_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
