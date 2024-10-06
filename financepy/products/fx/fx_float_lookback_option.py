##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.math import N
from ...utils.global_vars import g_days_in_year, g_small
from ...utils.error import FinError
from ...models.gbm_process_simulator import get_paths_times
from ...products.fx.fx_option import FXOption
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
# FLOAT STRIKE LOOKBACK CALL PAYS MAX(S(T)-s_min,0)
# FLOAT STRIKE LOOKBACK PUT PAYS MAX(s_max-S(T),0)
##########################################################################


class FXFloatLookbackOption(FXOption):
    """This is an FX option in which the strike of the option is not fixed
    but is set at expiry to equal the minimum fx rate in the case of a call
    or the maximum fx rate in the case of a put."""

    def __init__(self, expiry_dt: Date, option_type: OptionTypes):
        """Create the FX Float Look Back Option by specifying the expiry
        date and the option type."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.option_type = option_type

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
        """Valuation of the Floating Lookback option using Black-Scholes
        using the formulae derived by Goldman, Sosin and Gatto (1979)."""

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

        df = domestic_curve.df_t(t)
        r = -np.log(df) / t

        dq = foreign_curve.df_t(t)
        q = -np.log(dq) / t

        v = volatility
        s0 = stock_price
        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_min = stock_min_max
            if s_min > s0:
                raise FinError(
                    "s_min must be less than or equal to the stock price."
                )
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_max = stock_min_max
            if s_max < s0:
                raise FinError(
                    "s_max must be greater than or equal to the stock price."
                )

        if abs(r - q) < g_small:
            q = r + g_small

        dq = np.exp(-q * t)
        df = np.exp(-r * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = np.exp(b * t)

        # Taken from Haug Page 142
        if self.option_type == OptionTypes.EUROPEAN_CALL:

            a1 = (np.log(s0 / s_min) + (b + (v**2) / 2.0) * t) / v / np.sqrt(t)
            a2 = a1 - v * np.sqrt(t)

            if s_min == s0:
                term = N(-a1 + 2.0 * b * np.sqrt(t) / v) - expbt * N(-a1)
            elif s0 < s_min and w < -100:
                term = -expbt * N(-a1)
            else:
                term = ((s0 / s_min) ** (-w)) * N(
                    -a1 + 2.0 * b * np.sqrt(t) / v
                ) - expbt * N(-a1)

            v = s0 * dq * N(a1) - s_min * df * N(a2) + s0 * df * u * term

        elif self.option_type == OptionTypes.EUROPEAN_PUT:

            b1 = (np.log(s0 / s_max) + (b + (v**2) / 2.0) * t) / v / np.sqrt(t)
            b2 = b1 - v * np.sqrt(t)

            if s_max == s0:
                term = -N(b1 - 2.0 * b * np.sqrt(t) / v) + expbt * N(b1)
            elif s0 < s_max and w > 100:
                term = expbt * N(b1)
            else:
                term = (-((s0 / s_max) ** (-w))) * N(
                    b1 - 2.0 * b * np.sqrt(t) / v
                ) + expbt * N(b1)

            v = s_max * df * N(-b2) - s0 * dq * N(-b1) + s0 * df * u * term

        else:
            raise FinError(
                "Unknown lookback option type:" + str(self.option_type)
            )

        return v

    ##########################################################################

    def value_mc(
        self,
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min_max,
        num_paths=10000,
        num_steps_per_year=252,
        seed=4242,
    ):
        """Value FX floating lookback option using Monte Carlo"""
        t = (self.expiry_dt - value_dt) / g_days_in_year
        df = domestic_curve.df_t(t)
        r = -np.log(df) / t

        dq = foreign_curve.df_t(t)
        q = -np.log(dq) / t

        num_time_steps = int(t * num_steps_per_year)
        mu = r - q

        option_type = self.option_type
        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_min = stock_min_max
            if s_min > stock_price:
                raise FinError(
                    "s_min must be less than or equal to the stock price."
                )
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_max = stock_min_max
            if s_max < stock_price:
                raise FinError(
                    "s_max must be greater than or equal to the stock price."
                )

        t_all, s_all = get_paths_times(
            num_paths, num_time_steps, t, mu, stock_price, volatility, seed
        )

        # Due to anti-thetics we have doubled the number of paths
        num_paths = 2 * num_paths
        payoff = np.zeros(num_paths)

        if option_type == OptionTypes.EUROPEAN_CALL:
            # minimum on time dimension axis 1
            s_min_vector = np.min(s_all, axis=1)
            s_min_vector = np.minimum(s_min_vector, s_min)
            payoff = np.maximum(s_all[:, -1] - s_min_vector, 0.0)
        elif option_type == OptionTypes.EUROPEAN_PUT:
            # maximum on time dimension axis 1
            s_max_vector = np.max(s_all, axis=1)
            s_max_vector = np.maximum(s_max_vector, s_max)
            payoff = np.maximum(s_max_vector - s_all[:, -1], 0.0)
        else:
            raise FinError("Unknown lookback option type:" + str(option_type))

        v = payoff.mean() * np.exp(-r * t)
        return v


##########################################################################
