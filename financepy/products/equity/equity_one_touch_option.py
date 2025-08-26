##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np
from numba import njit

from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.global_types import TouchOptionTypes
from ...utils.error import FinError
from ...products.equity.equity_option import EquityOption
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve
from ...models.gbm_process_simulator import get_paths_times

from ...utils.math import normcdf_vect

########################################################################################
# TODO: Implement Sobol random numbers
# TODO: Improve convergence
########################################################################################

DEBUG_MODE = False


@njit(fastmath=True, cache=True)
def _barrier_pay_one_at_hit_pv_down(s, hh, r, dt):
    """Pay $1 if the stock crosses the barrier hh from above. PV payment."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] <= hh:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_one_at_hit_pv_up(s, hh, r, dt):
    """Pay $1 if the stock crosses the barrier hh from below. PV payment."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] >= hh:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_asset_at_expiry_down_out(s, hh):
    """Pay $1 if the stock crosses the barrier hh from above. PV payment."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] <= hh:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_asset_at_expiry_up_out(s, hh):
    """Pay $1 if the stock crosses the barrier hh from below. PV payment."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in range(0, num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] >= hh:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


########################################################################################


class EquityOneTouchOption(EquityOption):
    """A EquityOneTouchOption is an option in which the buyer receives one
    unit of cash OR stock if the stock price touches a barrier at any time
    before the option expiry date and zero otherwise. The choice of cash or
    stock is made at trade initiation. The single barrier payoff must define
    whether the option pays or cancels if the barrier is touched and also when
    the payment is made (at hit time or option expiry). All of these variants
    are all members of the FinTouchOptionTypes enumerated type."""

    def __init__(
        self,
        expiry_dt: Date,
        opt_type: TouchOptionTypes,
        barrier_price: float,
        payment_size: float = 1.0,
    ):
        """Create the one touch option by defining its expiry date and the
        barrier level and a payment size if it is a cash ."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.opt_type = opt_type
        self.barrier_price = float(barrier_price)
        self.payment_size = payment_size

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Equity One-Touch Option valuation using the Black-Scholes model
        assuming a continuous (American) barrier from value date to expiry.
        Handles both cash-or-nothing and asset-or-nothing options."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError("Discount Curve date not same as option valuation date")

        if dividend_curve.value_dt != value_dt:
            raise FinError("Dividend Curve date not same as option valuation date")

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        t = max(t, 1e-6)

        s0 = stock_price
        hh = self.barrier_price
        k = self.payment_size

        sqrt_t = np.sqrt(t)

        df = discount_curve.df(self.expiry_dt)
        r = discount_curve.cc_rate(self.expiry_dt)
        q = dividend_curve.cc_rate(self.expiry_dt)

        v = model.volatility
        v = max(v, 1e-6)

        # Using notation in Haug page 177
        b = r - q
        mu = (b - v * v / 2.0) / v / v
        lam = np.sqrt(mu * mu + 2.0 * r / v / v)

        if DEBUG_MODE:
            print("t:", t)
            print("vol", v)
            print("b", b)
            print("mu", mu)
            print("lam", lam)

        # Reference Option Pricing Formulas by Espen Gaarder Haug. Page 176.

        if self.opt_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT:
            # HAUG 1

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = 1.0
            z = np.log(hh / s0) / v / sqrt_t + lam * v * sqrt_t
            a5_1 = np.power(hh / s0, mu + lam) * normcdf_vect(eta * z)
            a5_2 = np.power(hh / s0, mu - lam) * normcdf_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (a5_1 + a5_2) * k
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_CASH_AT_HIT:
            # HAUG 2

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            z = np.log(hh / s0) / v / sqrt_t + lam * v * sqrt_t
            a5_1 = np.power(hh / s0, mu + lam) * normcdf_vect(eta * z)
            a5_2 = np.power(hh / s0, mu - lam) * normcdf_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (a5_1 + a5_2) * k
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT:
            # HAUG 3

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = 1.0
            k = hh
            z = np.log(hh / s0) / v / sqrt_t + lam * v * sqrt_t
            a5_1 = np.power(hh / s0, mu + lam) * normcdf_vect(eta * z)
            a5_2 = np.power(hh / s0, mu - lam) * normcdf_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (a5_1 + a5_2) * k
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT:
            # HAUG 4

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            k = hh
            z = np.log(hh / s0) / v / sqrt_t + lam * v * sqrt_t
            a5_1 = np.power(hh / s0, mu + lam) * normcdf_vect(eta * z)
            a5_2 = np.power(hh / s0, mu - lam) * normcdf_vect(
                eta * z - 2.0 * eta * lam * v * sqrt_t
            )
            v = (a5_1 + a5_2) * k
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY:
            # HAUG 5

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = +1.0
            phi = -1.0
            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b2 = k * df * normcdf_vect(phi * x2 - phi * v * sqrt_t)
            b4 = (
                k
                * df
                * np.power(hh / s0, 2.0 * mu)
                * normcdf_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b2 + b4
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY:
            # HAUG 6

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            phi = +1.0

            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b2 = k * df * normcdf_vect(phi * x2 - phi * v * sqrt_t)
            b4 = (
                k
                * df
                * np.power(hh / s0, 2.0 * mu)
                * normcdf_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b2 + b4
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 7

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = +1.0
            phi = -1.0
            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-q * t)
            a2 = s0 * dq * normcdf_vect(phi * x2)
            a4 = s0 * dq * np.power(hh / s0, 2.0 * (mu + 1.0)) * normcdf_vect(eta * y2)
            v = a2 + a4
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 8

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            phi = +1.0
            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-q * t)
            a2 = s0 * dq * normcdf_vect(phi * x2)
            a4 = s0 * dq * np.power(hh / s0, 2.0 * (mu + 1.0)) * normcdf_vect(eta * y2)
            v = a2 + a4
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING:
            # HAUG 9

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = +1.0
            phi = +1.0

            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b2 = k * df * normcdf_vect(phi * x2 - phi * v * sqrt_t)
            b4 = (
                k
                * df
                * np.power(hh / s0, 2.0 * mu)
                * normcdf_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b2 - b4
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING:
            # HAUG 10

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            phi = -1.0

            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            b2 = k * df * normcdf_vect(phi * x2 - phi * v * sqrt_t)
            b4 = (
                k
                * df
                * np.power(hh / s0, 2.0 * mu)
                * normcdf_vect(eta * y2 - eta * v * sqrt_t)
            )
            v = b2 - b4
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 11

            if np.any(s0 <= hh):
                raise FinError("Stock price is currently below barrier.")

            eta = +1.0
            phi = +1.0

            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-q * t)
            a2 = s0 * dq * normcdf_vect(phi * x2)
            a4 = s0 * dq * np.power(hh / s0, 2.0 * (mu + 1.0)) * normcdf_vect(eta * y2)
            v = a2 - a4
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 12

            if np.any(s0 >= hh):
                raise FinError("Stock price is currently above barrier.")

            eta = -1.0
            phi = -1.0

            x2 = np.log(s0 / hh) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            y2 = np.log(hh / s0) / v / sqrt_t + (mu + 1.0) * v * sqrt_t
            dq = np.exp(-q * t)
            a2 = s0 * dq * normcdf_vect(phi * x2)
            a4 = s0 * dq * np.power(hh / s0, 2.0 * (mu + 1.0)) * normcdf_vect(eta * y2)
            v = a2 - a4
            return v

        else:
            raise FinError("Unknown option type.")

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Touch Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Accuracy is not great when compared to the analytical
        result as we only observe the barrier a finite number of times. The
        convergence is slow."""

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t

        dq = dividend_curve.df(self.expiry_dt)
        q = -np.log(dq) / t

        num_time_steps = int(t * num_steps_per_year) + 1
        dt = t / num_time_steps

        v = model.volatility
        s0 = stock_price
        mu = r - q

        time_grid, s = get_paths_times(num_paths, num_time_steps, t, mu, s0, v, seed)

        hh = self.barrier_price
        xx = self.payment_size

        v = 0.0

        if self.opt_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT:
            # HAUG 1

            if s0 <= hh:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_down(s, hh, r, dt)
            v = v * xx
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_CASH_AT_HIT:
            # HAUG 2

            if s0 >= hh:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_up(s, hh, r, dt)
            v = v * xx
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT:
            # HAUG 3

            if s0 <= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_down(s, hh, r, dt) * hh
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT:
            # HAUG 4

            if s0 >= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_up(s, hh, r, dt) * hh
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY:
            # HAUG 5

            if s0 <= hh:
                raise FinError("Barrier has  ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_down(s, hh, 0.0, dt)
            v = v * xx * np.exp(-r * t)
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY:
            # HAUG 6

            if s0 >= hh:
                raise FinError("Barrier has ALREADY been crossed.")

            v = _barrier_pay_one_at_hit_pv_up(s, hh, 0.0, dt)
            v = v * xx * np.exp(-r * t)
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 7

            if s0 <= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_down(s, hh, 0.0, dt) * hh
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY:
            # HAUG 8

            if s0 >= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_one_at_hit_pv_up(s, hh, 0.0, dt) * hh
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING:
            # HAUG 9

            if s0 <= hh:
                raise FinError("Barrier has ALREADY been crossed.")

            v = 1.0 - _barrier_pay_one_at_hit_pv_down(s, hh, 0.0, dt)
            v = v * xx * np.exp(-r * t)
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING:
            # HAUG 10

            if s0 >= hh:
                raise FinError("Barrier has ALREADY been crossed.")

            v = 1.0 - _barrier_pay_one_at_hit_pv_up(s, hh, 0.0, dt)
            v = v * xx * np.exp(-r * t)
            return v

        elif self.opt_type == TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 11

            if s0 <= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_asset_at_expiry_down_out(s, hh)
            v = v * np.exp(-r * t)
            return v

        elif self.opt_type == TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING:
            # HAUG 12

            if s0 >= hh:
                raise FinError("Stock price is currently below barrier.")

            v = _barrier_pay_asset_at_expiry_up_out(s, hh)
            v = v * np.exp(-r * t)
            return v
        else:
            raise FinError("Unknown option type.")

        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("OPTION TYPE", self.opt_type)
        s += label_to_string("BARRIER LEVEL", self.barrier_price)
        s += label_to_string("PAYMENT SIZE", self.payment_size, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
