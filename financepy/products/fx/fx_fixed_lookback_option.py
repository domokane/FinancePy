##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...utils.math import N
from ...utils.global_vars import g_days_in_year, g_small
from ...utils.error import FinError
from ...models.gbm_process_simulator import get_paths_times
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.global_types import OptionTypes
from ...market.curves.discount_curve import DiscountCurve

##########################################################################
# TODO: Attempt control variate adjustment to monte carlo
# TODO: Sobol for Monte Carlo
# TODO: TIGHTEN UP LIMIT FOR W FROM 100
# TODO: Vectorise the analytical pricing formula
##########################################################################

##########################################################################
# FIXED STRIKE LOOKBACK CALL PAYS MAX(SMAX-K,0)
# FIXED STRIKE LOOKBACK PUT PAYS MAX(K-SMIN,0)
##########################################################################


class FXFixedLookbackOption:
    """The Class for FX Fixed Strike Lookback options."""

    def __init__(
        self, expiry_dt: Date, option_type: OptionTypes, option_strike: float
    ):
        """Create option with expiry date, option type and the option strike"""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.option_type = option_type
        self.option_strike = option_strike

    ##########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        volatility: float,
        stock_min_max: float,
    ):
        """Value FX Fixed Lookback Option using Black Scholes model and
        analytical formulae."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as option value date"
            )

        if foreign_curve.value_dt != value_dt:
            raise FinError(
                "Foreign Curve valuation date not same as option value date"
            )

        t = (self.expiry_dt - value_dt) / g_days_in_year

        df = domestic_curve.df(self.expiry_dt)
        r = -np.log(df) / t

        dq = foreign_curve.df(self.expiry_dt)
        q = -np.log(dq) / t

        v = volatility
        s0 = stock_price
        k = self.option_strike
        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_max = stock_min_max
            if s_max < s0:
                raise FinError("The Smax value must be >= the stock price.")
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_min = stock_min_max
            if s_min > s0:
                raise FinError("The s_min value must be <= the stock price.")

        # There is a risk of an overflow in the limit of q=r which
        # we remove by adjusting the value of the dividend
        if abs(r - q) < g_small:
            q = r + g_small

        df = exp(-r * t)
        dq = exp(-q * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = exp(b * t)

        # Taken from Hull Page 536 (6th edition) and Haug Page 143
        if self.option_type == OptionTypes.EUROPEAN_CALL:

            if k > s_max:
                d1 = (log(s0 / k) + (b + v * v / 2.0) * t) / v / sqrt(t)
                d2 = d1 - v * sqrt(t)

                if s0 == k:
                    term = -N(d1 - 2.0 * b * sqrt(t) / v) + expbt * N(d1)
                elif s0 < k and w > 100:
                    term = expbt * N(d1)
                else:
                    term = (-((s0 / k) ** (-w))) * N(
                        d1 - 2.0 * b * sqrt(t) / v
                    ) + expbt * N(d1)

                v = s0 * dq * N(d1) - k * df * N(d2) + s0 * df * u * term

            else:
                e1 = (log(s0 / s_max) + (b + v * v / 2.0) * t) / v / sqrt(t)
                e2 = e1 - v * sqrt(t)

                if s0 == s_max:
                    term = -N(e1 - 2.0 * b * sqrt(t) / v) + expbt * N(e1)
                elif s0 < s_max and w > 100:
                    term = expbt * N(e1)
                else:
                    term = (-((s0 / s_max) ** (-w))) * N(
                        e1 - 2.0 * b * sqrt(t) / v
                    ) + expbt * N(e1)

                v = (
                    df * (s_max - k)
                    + s0 * dq * N(e1)
                    - s_max * df * N(e2)
                    + s0 * df * u * term
                )

        elif self.option_type == OptionTypes.EUROPEAN_PUT:

            if k >= s_min:
                f1 = (log(s0 / s_min) + (b + v * v / 2.0) * t) / v / sqrt(t)
                f2 = f1 - v * sqrt(t)

                if s0 == s_min:
                    term = N(-f1 + 2.0 * b * sqrt(t) / v) - expbt * N(-f1)
                elif s0 > s_min and w < -100:
                    term = -expbt * N(-f1)
                else:
                    term = ((s0 / s_min) ** (-w)) * N(
                        -f1 + 2.0 * b * sqrt(t) / v
                    ) - expbt * N(-f1)

                v = (
                    df * (k - s_min)
                    - s0 * dq * N(-f1)
                    + s_min * df * N(-f2)
                    + s0 * df * u * term
                )

            else:
                d1 = (log(s0 / k) + (b + v * v / 2) * t) / v / sqrt(t)
                d2 = d1 - v * sqrt(t)
                if s0 == k:
                    term = N(-d1 + 2.0 * b * sqrt(t) / v) - expbt * N(-d1)
                elif s0 > k and w < -100:
                    term = -expbt * N(-d1)
                else:
                    term = ((s0 / k) ** (-w)) * N(
                        -d1 + 2.0 * b * sqrt(t) / v
                    ) - expbt * N(-d1)

                v = k * df * N(-d2) - s0 * dq * N(-d1) + s0 * df * u * term

        else:
            raise FinError(
                "Unknown lookback option type:" + str(self.option_type)
            )

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float,  # FORDOM
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        volatility: float,
        spot_fx_rate_min_max: float,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Value FX Fixed Lookback option using Monte Carlo."""

        t = (self.expiry_dt - value_dt) / g_days_in_year
        s_0 = spot_fx_rate

        df = domestic_curve.df_t(t)
        r_d = -np.log(df) / t

        dq = foreign_curve.df_t(t)
        r_f = -np.log(dq) / t

        mu = r_d - r_f

        num_time_steps = int(t * num_steps_per_year)

        option_type = self.option_type
        k = self.option_strike

        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_max = spot_fx_rate_min_max
            if s_max < s_0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price."
                )
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_min = spot_fx_rate_min_max
            if s_min > s_0:
                raise FinError(
                    "s_min must be less than or equal to the stock price."
                )

        t_all, s_all = get_paths_times(
            num_paths, num_time_steps, t, mu, s_0, volatility, seed
        )

        payoff = np.zeros(num_paths)

        if option_type == OptionTypes.EUROPEAN_CALL:
            s_max_vector = np.max(s_all, axis=1)
            s_maxs = np.ones(num_paths) * s_max
            payoff = np.maximum(s_max_vector - k, 0.0)
            payoff = np.maximum(payoff, s_maxs - k)
        elif option_type == OptionTypes.EUROPEAN_PUT:
            s_min_vector = np.min(s_all, axis=1)
            s_mins = np.ones(num_paths) * s_min
            payoff = np.maximum(k - s_min_vector, 0.0)
            payoff = np.maximum(payoff, k - s_mins)
        else:
            raise FinError("Unknown lookback option type:" + str(option_type))

        v = payoff.mean() * exp(-r_d * t)
        return v
