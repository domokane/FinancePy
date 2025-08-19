##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np


from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.error import FinError
from ...utils.global_types import TouchOptionTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve
from ...models.gbm_process_simulator import get_paths_times
from ...products.fx.fx_option import FXOption

from numba import njit

from ...utils.math import n_vect

###############################################################################
# TODO: Implement Sobol random numbers
# TODO: Improve convergence
# TODO: Fix risk numbers
###############################################################################

###############################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_one_at_hit_pv_down(s, H, r, dt):
    """Pay $1 if the stock crosses the barrier H from above. PV payment."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] <= H:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


###############################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_one_at_hit_pv_up(s, H, r, dt):
    """Pay $1 if the stock crosses the barrier H from below. PV payment."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] >= H:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


###############################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_asset_at_expiry_down_out(s, h):
    """Pay $1 if the stock crosses the barrier H from above. PV payment."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] <= h:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


###############################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_asset_at_expiry_up_out(s, h):
    """Pay $1 if the stock crosses the barrier H from below. PV payment."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] >= h:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


###############################################################################


class FXOneTouchOption(FXOption):
    """A FinFXOneTouchOption is an option in which the buyer receives one
    unit of currency if the FX rate touches a barrier at any time
    before the option expiry date and zero otherwise. The single barrier
    payoff must define whether the option pays or cancels if the barrier is
    touched and also when the payment is made (at hit time or option expiry).
    All of these variants are members of the FinTouchOptionTypes type."""

    def __init__(
        self,
        expiry_dt: Date,
        option_type: TouchOptionTypes,
        barrier_rate: float,
        payment_size: float = 1.0,
    ):
        """Create the one touch option by defining its expiry date and the
        barrier level and a payment size if it is a cash ."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.option_type = option_type
        self.barrier_rate = float(barrier_rate)
        self.payment_size = payment_size

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: Union[float, np.ndarray],
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        model,
    ):
        """FX One-Touch Option valuation using the Black-Scholes model
        assuming a continuous (American) barrier from value date to expiry.
        Handles both cash-or-nothing and asset-or-nothing options."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError("Domestic Curve date not same as valuation date")

        if foreign_curve.value_dt != value_dt:
            raise FinError("Foreign Curve date not same as valuation date")

        debug_mode = False

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        t = max(t, 1e-6)

        s0 = spot_fx_rate
        h = self.barrier_rate
        k = self.payment_size

        sqrt_t = np.sqrt(t)

        df = domestic_curve.df(self.expiry_dt)
        r_d = domestic_curve.cc_rate(self.expiry_dt)
        r_f = foreign_curve.cc_rate(self.expiry_dt)

        v = model.volatility
        v = max(v, 1e-6)

        # Using notation in Haug page 177
        b = r_d - r_f
        mu = (b - v * v / 2.0) / v / v
        lam = np.sqrt(mu * mu + 2.0 * r_d / v / v)

        if debug_mode:
            print("t:", t)
            print("vol", v)
            print("b", b)
            print("mu", mu)
            print("lam", lam)

        if self.option_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT:
            # HAUG 1

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = 1.0
            z = np.log(h / s0) / v / sqrt_t + lam * v * sqrt_t
            A5_1 = np.power(h / s0, mu + lam) * n_vect(eta * z)
            A5_2 = np.power(h / s0, mu - lam) * n_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (A5_1 + A5_2) * k
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_CASH_AT_HIT:
            # HAUG 2

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            z = np.log(h / s0) / v / sqrt_t + lam * v * sqrt_t
            A5_1 = np.power(h / s0, mu + lam) * n_vect(eta * z)
            A5_2 = np.power(h / s0, mu - lam) * n_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (A5_1 + A5_2) * k
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT:
            # HAUG 3

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = 1.0
            k = h
            z = np.log(h / s0) / v / sqrt_t + lam * v * sqrt_t
            A5_1 = np.power(h / s0, mu + lam) * n_vect(eta * z)
            A5_2 = np.power(h / s0, mu - lam) * n_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (A5_1 + A5_2) * k
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT:
            # HAUG 4

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            k = h
            z = np.log(h / s0) / v / sqrt_t + lam * v * sqrt_t
            A5_1 = np.power(h / s0, mu + lam) * n_vect(eta * z)
            A5_2 = np.power(h / s0, mu - lam) * n_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (A5_1 + A5_2) * k
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY:
            # HAUG 5

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = +1.0
            phi = -1.0
            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b_2 = k * df * n_vect(phi * x2 - phi * v * sqrt_t)
            b_4 = (
                k
                * df
                * np.power(h / s0, 2.0 * mu)
                * n_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b_2 + b_4
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY:
            # HAUG 6

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            phi = +1.0

            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b_2 = k * df * n_vect(phi * x2 - phi * v * sqrt_t)
            b_4 = (
                k
                * df
                * np.power(h / s0, 2.0 * mu)
                * n_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b_2 + b_4
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 7

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = +1.0
            phi = -1.0
            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-r_f * t)
            a_2 = s0 * dq * n_vect(phi * x2)
            a_4 = s0 * dq * np.power(h / s0, 2.0 * (mu + 1.0)) * n_vect(eta * y2)
            v = a_2 + a_4
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 8

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            phi = +1.0
            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-r_f * t)
            a_2 = s0 * dq * n_vect(phi * x2)
            a_4 = s0 * dq * np.power(h / s0, 2.0 * (mu + 1.0)) * n_vect(eta * y2)
            v = a_2 + a_4
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING:
            # HAUG 9

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = +1.0
            phi = +1.0

            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b_2 = k * df * n_vect(phi * x2 - phi * v * sqrt_t)
            b_4 = (
                k
                * df
                * np.power(h / s0, 2.0 * mu)
                * n_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b_2 - b_4
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING:
            # HAUG 10

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            phi = -1.0

            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b_2 = k * df * n_vect(phi * x2 - phi * v * sqrt_t)
            b_4 = (
                k
                * df
                * np.power(h / s0, 2.0 * mu)
                * n_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b_2 - b_4
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 11

            if np.any(s0 <= h):
                raise FinError("FX Rate is currently below barrier.")

            eta = +1.0
            phi = +1.0

            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-r_f * t)
            a_2 = s0 * dq * n_vect(phi * x2)
            a_4 = s0 * dq * np.power(h / s0, 2.0 * (mu + 1.0)) * n_vect(eta * y2)
            v = a_2 - a_4
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 12

            if np.any(s0 >= h):
                raise FinError("FX Rate is currently above barrier.")

            eta = -1.0
            phi = -1.0

            x2 = np.log(s0 / h) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(h / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-r_f * t)
            a_2 = s0 * dq * n_vect(phi * x2)
            a_4 = s0 * dq * np.power(h / s0, 2.0 * (mu + 1.0)) * n_vect(eta * y2)
            v = a_2 - a_4
            return v

        else:
            raise FinError("Unknown option type.")

        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("OPTION TYPE", self.option_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_rate)
        s += label_to_string("PAYMENT SIZE", self.payment_size, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        dom_curve: DiscountCurve,
        for_curve: DiscountCurve,
        model,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Touch Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Accuracy is not great when compared to the analytical
        result as we only observe the barrier a finite number of times. The
        convergence is slow."""

        # "THIS NEEDS TO BE CHECKED"

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        df_d = dom_curve.df(self.expiry_dt)
        r_d = -np.log(df_d) / t

        df_f = for_curve.df(self.expiry_dt)
        r_f = -np.log(df_f) / t

        num_time_steps = int(t * num_steps_per_year) + 1
        dt = t / num_time_steps

        v = model.volatility
        s0 = stock_price
        mu = r_d - r_f

        tgrid, s = get_paths_times(num_paths, num_time_steps, t, mu, s0, v, seed)

        h = self.barrier_rate
        x = self.payment_size

        v = 0.0

        if self.option_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT:
            # HAUG 1

            if s0 <= h:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_down(s, h, r_d, dt)
            v = v * x
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_CASH_AT_HIT:
            # HAUG 2

            if s0 >= h:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_up(s, h, r_d, dt)
            v = v * x
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT:
            # HAUG 3

            if s0 <= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_down(s, h, r_d, dt) * h
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT:
            # HAUG 4

            if s0 >= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_up(s, h, r_d, dt) * h
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY:
            # HAUG 5

            if s0 <= h:
                raise FinError("Barrier has  ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_down(s, h, 0.0, dt)
            v = v * x * np.exp(-r_d * t)
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY:
            # HAUG 6

            if s0 >= h:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_up(s, h, 0.0, dt)
            v = v * x * np.exp(-r_d * t)
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 7

            if s0 <= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_down(s, h, 0.0, dt) * h
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 8

            if s0 >= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_up(s, h, 0.0, dt) * h
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING:
            # HAUG 9

            if s0 <= h:
                raise FinError("Barrier has ALREADY been crossed.")

            v = 1.0 - _barrier_pay_one_at_hit_pv_down(s, h, 0.0, dt)
            v = v * x * np.exp(-r_d * t)
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING:
            # HAUG 10

            if s0 >= h:
                raise FinError("Barrier has ALREADY been crossed.")

            v = 1.0 - _barrier_pay_one_at_hit_pv_up(s, h, 0.0, dt)
            v = v * x * np.exp(-r_d * t)
            return v

        elif self.option_type == TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 11

            if s0 <= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_asset_at_expiry_down_out(s, h)
            v = v * np.exp(-r_d * t)
            return v

        elif self.option_type == TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 12

            if s0 >= h:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_asset_at_expiry_up_out(s, h)
            v = v * np.exp(-r_d * t)
            return v
        else:
            raise FinError("Unknown option type.")

        return v


###############################################################################
