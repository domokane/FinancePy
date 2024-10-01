##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...utils.math import N
from ...utils.global_vars import g_days_in_year, g_small
from ...utils.error import FinError
from ...utils.date import Date

from ...models.gbm_process_simulator import get_paths_times
from ...products.equity.equity_option import EquityOption
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from ...utils.global_types import OptionTypes

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


class EquityFixedLookbackOption(EquityOption):
    """This is an equity option in which the strike of the option is fixed but
    the value of the stock price used to determine the payoff is the maximum
    in the case of a call option, and a minimum in the case of a put option."""

    def __init__(
        self, expiry_dt: Date, option_type: OptionTypes, strike_price: float
    ):
        """Create the FixedLookbackOption by specifying the expiry date, the
        option type and the option strike."""

        check_argument_types(self.__init__, locals())

        if (
            option_type != OptionTypes.EUROPEAN_CALL
            and option_type != OptionTypes.EUROPEAN_PUT
        ):
            raise FinError("Option type must be EUROPEAN_CALL or EUROPEAN_PUT")

        self.expiry_dt = expiry_dt
        self.option_type = option_type
        self.strike_price = strike_price

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        volatility: float,
        stock_min_max: float,
    ):
        """Valuation of the Fixed Lookback option using Black-Scholes using
        the formulae derived by Conze and Viswanathan (1991). One of the inputs
        is the minimum of maximum of the stock price since the start of the
        option depending on whether the option is a call or a put."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if dividend_curve.value_dt != value_dt:
            raise FinError(
                "Dividend Curve valuation date not same as option value date"
            )

        t = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t

        v = volatility
        s0 = stock_price
        k = self.strike_price
        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_max = stock_min_max
            if s_max < s0:
                raise FinError("The Smax value must be >= the stock price.")
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_min = stock_min_max
            if s_min > s0:
                raise FinError("The Smin value must be <= the stock price.")

        # There is a risk of an overflow in the limit of q=r which
        # we remove by adjusting the value of the dividend
        if abs(r - q) < g_small:
            q = r + g_small

        df = np.exp(-r * t)
        dq = np.exp(-q * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / (v * v)
        expbt = np.exp(b * t)
        sqrt_t = np.sqrt(t)

        # Taken from Hull Page 536 (6th edition) and Haug Page 143
        if self.option_type == OptionTypes.EUROPEAN_CALL:

            if k > s_max:

                d1 = (np.log(s0 / k) + (b + v * v / 2.0) * t) / v / sqrt_t
                d2 = d1 - v * sqrt_t

                if s0 == k:
                    term = -N(d1 - 2.0 * b * sqrt_t / v) + expbt * N(d1)
                elif s0 < k and w > 100.0:
                    term = expbt * N(d1)
                else:
                    term = -np.power(s0 / k, -w) * N(
                        d1 - 2 * b * sqrt_t / v
                    ) + expbt * N(d1)

                v = s0 * dq * N(d1) - k * df * N(d2) + s0 * df * u * term

            else:

                e1 = (
                    (np.log(s0 / s_max) + (r - q + v * v / 2) * t) / v / sqrt_t
                )
                e2 = e1 - v * sqrt_t

                if s0 == s_max:
                    term = -N(e1 - 2.0 * b * sqrt_t / v) + expbt * N(e1)
                elif s0 < s_max and w > 100.0:
                    term = expbt * N(e1)
                else:
                    term = (-((s0 / s_max) ** (-w))) * N(
                        e1 - 2.0 * b * sqrt_t / v
                    ) + expbt * N(e1)

                v = (
                    df * (s_max - k)
                    + s0 * dq * N(e1)
                    - s_max * df * N(e2)
                    + s0 * df * u * term
                )

        elif self.option_type == OptionTypes.EUROPEAN_PUT:

            if k >= s_min:

                f1 = (np.log(s0 / s_min) + (b + v * v / 2.0) * t) / v / sqrt_t
                f2 = f1 - v * sqrt_t

                if s0 == s_min:
                    term = N(-f1 + 2.0 * b * sqrt_t / v) - expbt * N(-f1)
                elif s0 > s_min and w < -100.0:
                    term = -expbt * N(-f1)
                else:
                    term = ((s0 / s_min) ** (-w)) * N(
                        -f1 + 2.0 * b * sqrt_t / v
                    ) - expbt * N(-f1)

                v = (
                    df * (k - s_min)
                    - s0 * dq * N(-f1)
                    + s_min * df * N(-f2)
                    + s0 * df * u * term
                )

            else:

                d1 = (np.log(s0 / k) + (b + v * v / 2) * t) / v / sqrt_t
                d2 = d1 - v * sqrt_t

                if s0 == k:
                    term = N(-d1 + 2.0 * b * sqrt_t / v) - expbt * N(-d1)
                elif s0 > k and w < -100.0:
                    term = -expbt * N(-d1)
                else:
                    term = ((s0 / k) ** (-w)) * N(
                        -d1 + 2.0 * b * sqrt_t / v
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
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        volatility: float,
        stock_min_max: float,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Monte Carlo valuation of a fixed strike lookback option using a
        Black-Scholes model that assumes the stock follows a GBM process."""

        t = (self.expiry_dt - value_dt) / g_days_in_year

        df = discount_curve.df(self.expiry_dt)
        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        mu = r - q
        num_time_steps = int(t * num_steps_per_year)

        option_type = self.option_type
        k = self.strike_price

        s_min = 0.0
        s_max = 0.0

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            s_max = stock_min_max
            if s_max < stock_price:
                raise FinError(
                    "Smax must be greater than or equal to the stock price."
                )
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            s_min = stock_min_max
            if s_min > stock_price:
                raise FinError(
                    "Smin must be less than or equal to the stock price."
                )

        t_all, s_all = get_paths_times(
            num_paths, num_time_steps, t, mu, stock_price, volatility, seed
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

        v = payoff.mean() * df
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.option_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
