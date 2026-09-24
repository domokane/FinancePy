########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from typing import Union
import numpy as np

from ...utils.global_types import DoubleBarrierTypes

from ...utils.global_vars import G_DAYS_IN_YEAR
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve
from ...models.gbm_process_simulator import get_paths_times
from .fx_option import FXOption
from ...utils.check_values import check_curve_dt
from ...utils.helpers import option_years
from ...models.double_touch_option import fast_double_no_touch_pricer
from ...models.double_touch_option import barrier_pay_at_expiry_double_hit
from ...models.double_touch_option import p_double_touch_bb_parallel

########################################################################################
# TODO: Implement Sobol random numbers
# TODO: Improve convergence
# TODO: Fix risk numbers
########################################################################################


class FXDoubleOneTouchOption(FXOption):
    """A FinFXOneTouchOption is an option in which the buyer receives a
    rebate if the FX rate touches either a lower barrier or an upper barrier
    at any time before the option expiry date and zero otherwise."""

    def __init__(
        self,
        expiry_dt: Date,
        option_type: DoubleBarrierTypes,
        lower_barrier_fx_rate: float,
        upper_barrier_fx_rate: float,
        payment_size: float = 1.0,
    ) -> None:
        """Create the double one touch option by defining its expiry date and the
        barrier level and a payment size if it is a cash ."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.option_type = option_type
        self.lower_barrier_fx_rate = float(lower_barrier_fx_rate)
        self.upper_barrier_fx_rate = float(upper_barrier_fx_rate)
        self.payment_size = payment_size

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: Union[float, np.ndarray],
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
    ):
        """FX One-Touch Option valuation using the Black-Scholes model
        assuming a continuous (American) barrier from value date to expiry.
        Handles both cash-or-nothing and asset-or-nothing options."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        df_d = domestic_curve.df(self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        if spot_fx_rate < self.lower_barrier_fx_rate:
            v = self.payment_size * df_d
            return v

        if spot_fx_rate > self.upper_barrier_fx_rate:
            v = self.payment_size * df_d
            return v

        opt_type_value = self.option_type.value

        # price following Haug 4.89
        S0 = float(spot_fx_rate)
        L = self.lower_barrier_fx_rate
        U = self.upper_barrier_fx_rate
        K = float(self.payment_size)
        sigma = float(model.volatility)

        t_exp = option_years(value_dt, self.expiry_dt)

        # moved code to Numba for speed
        c = fast_double_no_touch_pricer(S0, L, U, K, t_exp, opt_type_value, r_d, r_f, sigma)

        if self.option_type == DoubleBarrierTypes.KNOCK_OUT:
            pv = c
        else:
            pv = K * df_d - c

        return pv

    ###########################################################################

    def value_mc_slow(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Double one touch Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Accuracy is not great when compared to the analytical
        result as we only observe the barrier a finite number of times. The
        convergence is slow."""

        t_exp = option_years(value_dt, self.expiry_dt)

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        df_d = domestic_curve.df(self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        num_time_steps = int(t_exp * num_steps_per_year) + 1

        v = model.volatility
        s0 = stock_price
        mu = r_d - r_f

        _, s = get_paths_times(num_paths, num_time_steps, t_exp, mu, s0, v, seed)

        k1 = self.lower_barrier_fx_rate
        k2 = self.upper_barrier_fx_rate

        p_hit = barrier_pay_at_expiry_double_hit(s, k1, k2)

        if self.option_type == DoubleBarrierTypes.KNOCK_OUT:
            pv = self.payment_size * (1.0 - p_hit) * df_d
        else:
            pv = self.payment_size * p_hit * df_d

        return pv

    ####################################################################################

    def value_mc_bb_slow(
        self,
        value_dt: Date,
        stock_price: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_steps_per_year: int=52,
        num_paths: int=100000,
        seed=42,
    ):

        # Double one-touch option valuation with Brownian-bridge correction.
        t_exp = option_years(value_dt, self.expiry_dt)

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        df_d = domestic_curve.df(self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        mu = r_d - r_f
        sigma = model.volatility
        s0 = stock_price

        np.random.seed(seed)

        num_steps = int(num_steps_per_year * t_exp)
        dt = t_exp / num_steps
        nudt = (mu - 0.5 * sigma * sigma) * dt
        sigsdt = sigma * np.sqrt(dt)

        lnL = np.log(self.lower_barrier_fx_rate)
        lnU = np.log(self.upper_barrier_fx_rate)

        hits = 0

        for _ in range(num_paths):
            x = np.log(s0)
            hit = False
            for _ in range(num_steps):
                x_new = x + nudt + sigsdt * np.random.randn()

                # sure hits at endpoints
                if x <= lnL or x >= lnU or x_new <= lnL or x_new >= lnU:
                    hit = True
                    break

                # Brownian bridge correction
                # Lower barrier
                if x > lnL and x_new > lnL:
                    p_hit_lower = np.exp(-2.0 * (x - lnL) * (x_new - lnL) / (sigma**2 * dt))
                    if np.random.rand() < p_hit_lower:
                        hit = True
                        break

                # Upper barrier
                if x < lnU and x_new < lnU:
                    p_hit_upper = np.exp(-2.0 * (lnU - x) * (lnU - x_new) / (sigma**2 * dt))
                    if np.random.rand() < p_hit_upper:
                        hit = True
                        break

                x = x_new

            if hit:
                hits += 1

        p_touch = hits / num_paths

        if self.option_type == DoubleBarrierTypes.KNOCK_OUT:
            # Double no-touch (KO): pays if not touched
            pv = self.payment_size * df_d * (1.0 - p_touch)
        else:
            # Double one-touch (KI): pays if touched
            pv = self.payment_size * df_d * p_touch

        return pv

    ####################################################################################

    def value_mc(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model: Model,
        num_steps_per_year: int = 52,
        num_paths: int = 1_000_000,
        seed: int = 42,
    ) -> float:
        """
        PV via Brownian-bridge Monte Carlo with Numba parallelism.
        - KNOCK_IN  (double one-touch): pays K if a barrier is touched before expiry.
        - KNOCK_OUT (double no-touch):   pays K if neither barrier is touched before expiry.
        """
        T = (self.expiry_dt - value_dt) / G_DAYS_IN_YEAR
        if T <= 0.0:
            return 0.0

        check_curve_dt(value_dt, domestic_curve)
        check_curve_dt(value_dt, foreign_curve)

        df_d = domestic_curve.df(self.expiry_dt)
        r_d = domestic_curve.zero_rate_cc(self.expiry_dt)
        r_f = foreign_curve.zero_rate_cc(self.expiry_dt)

        mu = r_d - r_f
        sigma = float(model.volatility)
        S0 = spot_fx_rate
        L = self.lower_barrier_fx_rate
        U = self.upper_barrier_fx_rate

        steps = int(max(1, num_steps_per_year * T))

        # Estimate touch probability with Brownian bridge (parallel)
        p_touch = p_double_touch_bb_parallel(S0, L, U, mu, sigma, T, steps, num_paths, seed)

        if self.option_type == DoubleBarrierTypes.KNOCK_OUT:
            # Double no-touch
            p = 1.0 - p_touch
        else:
            # Double one-touch
            p = p_touch

        return self.payment_size * df_d * p

    ####################################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("LOWER BARRIER FX RATE", self.lower_barrier_fx_rate)
        s += label_to_string("UPPER BARRIER FX RATE", self.upper_barrier_fx_rate)
        s += label_to_string("PAYMENT SIZE", self.payment_size, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
