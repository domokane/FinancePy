########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from typing import Union
import numpy as np
from numba import njit, prange

from ...utils.global_types import DoubleBarrierTypes

from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.error import FinError
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve
from ...models.gbm_process_simulator import get_paths_times
from .fx_option import FXOption

########################################################################################
# TODO: Implement Sobol random numbers
# TODO: Improve convergence
# TODO: Fix risk numbers
########################################################################################


@njit(fastmath=True, parallel=True, cache=True)
def _p_touch_bb_parallel(
    S0: float,
    L: float,
    U: float,
    mu: float,
    sigma: float,
    T: float,
    steps: int,
    n_paths: int,
    seed: int,
) -> float:
    """
    Estimate P(touch before T) for a double barrier using Brownian bridge.
    Parallel across paths; equality to a barrier counts as a hit.
    """
    if steps < 1:
        steps = 1

    logS0 = np.log(S0)
    lnL = np.log(L)
    lnU = np.log(U)

    dt = T / steps
    nudt = (mu - 0.5 * sigma * sigma) * dt
    sigsdt = sigma * np.sqrt(dt)
    sig2_dt_inv = 1.0 / (sigma * sigma * dt)

    hits = 0

    # Each path uses its own seed offset for reproducibility
    for p in prange(n_paths):
        # Simple per-path reseeding; good enough in practice
        np.random.seed(seed + 1315423911 ^ (p * 2654435761))
        x = logS0
        hit = False

        # immediate hit at t=0?
        if x <= lnL or x >= lnU:
            hits += 1
            continue

        for _ in range(steps):
            z = np.random.randn()
            x_new = x + nudt + sigsdt * z

            # endpoint breaches (include equality as hit)
            if (x <= lnL) or (x >= lnU) or (x_new <= lnL) or (x_new >= lnU):
                hit = True
                break

            # Brownian-bridge hit probability between x and x_new
            # Lower barrier (a = lnL), valid when both endpoints strictly above lnL
            p_lower = 0.0
            if (x > lnL) and (x_new > lnL):
                t1 = (x - lnL) * (x_new - lnL)
                # guard tiny negatives due to rounding
                expo = -2.0 * t1 * sig2_dt_inv
                if expo < 0.0:
                    p_lower = np.exp(expo)
                else:
                    p_lower = 1.0  # if expo >= 0 due to numerical edge, force hit

            # Upper barrier (b = lnU), valid when both endpoints strictly below lnU
            p_upper = 0.0
            if (x < lnU) and (x_new < lnU):
                t2 = (lnU - x) * (lnU - x_new)
                expo = -2.0 * t2 * sig2_dt_inv
                if expo < 0.0:
                    p_upper = np.exp(expo)
                else:
                    p_upper = 1.0

            # Combined probability of hitting either barrier in (t, t+dt)
            # Assuming independence conditional on endpoints (standard approximation)
            p_any = 1.0 - (1.0 - p_lower) * (1.0 - p_upper)
            u = np.random.rand()
            if u < p_any:
                hit = True
                break

            x = x_new

        if hit:
            hits += 1

    return hits / n_paths


########################################################################################


@njit(fastmath=True, cache=True)
def _barrier_pay_at_expiry_double_hit(s, k1, k2):
    """Pay $1 if the stock crosses the barrier H from above. PV payment."""
    num_paths, num_time_steps = s.shape
    hits = 0.0

    for ip in range(0, num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            x = s[ip][it]
            if not (k1 < x < k2):
                hit_flag = 1
                break

        hits = hits + hit_flag

    p_hit = hits / num_paths
    return p_hit


########################################################################################


@njit(fastmath=True, cache=True)
def _fast_double_no_touch_pricer(S0, L, U, K, t_exp, opt_type, r_d, r_f, sigma):

    df_d = np.exp(-r_d * t_exp)

    # Immediate touch
    if S0 <= L or S0 >= U:
        if opt_type == 2:  # DNT
            return 0.0
        else:  # DOT
            return K * df_d

    # Precompute constants
    Z = np.log(U / L)
    b = r_d - r_f
    sig2 = sigma * sigma
    term1 = 2.0 * b / sig2 - 1.0
    alpha = -0.5 * term1
    beta = -0.25 * term1 * term1 - 2.0 * r_d / sig2

    logSL = np.log(S0 / L)
    logSU = np.log(S0 / U)

    # Heuristic n_max from damping bound
    # exp(-0.5*sig2*(n*pi/Z)^2 * t) <= eps  => n >= ...
    eps = 1e-14
    base = (Z / (np.pi * sigma * np.sqrt(2.0 * t_exp))) * np.sqrt(np.log(1.0 / eps))
    n_max = int(base) + 5
    if n_max < 50:
        n_max = 50
    if n_max > 2000:
        n_max = 2000

    # Series sum
    c = 0.0
    rel_tol = 1e-12  # finance-grade precision

    for i in range(1, n_max + 1):
        # m = i * pi / Z
        m = (i * np.pi) / Z

        # (-1)^i without pow
        alt = -1.0 if (i & 1) else 1.0

        # (S0/L)^alpha and (S0/U)^alpha
        SL_alpha = np.exp(alpha * logSL)
        SU_alpha = np.exp(alpha * logSU)

        denom = (alpha * alpha) + (m * m)
        term_i_b = SL_alpha - alt * SU_alpha

        # sin(m * log(S0/L))
        s = np.sin(m * logSL)

        # exp damp: exp(-0.5*(m^2 - beta)*sig2*t)
        damp = np.exp(-0.5 * (m * m - beta) * sig2 * t_exp)

        # assemble term: 2*pi*i*K/Z^2 * (term_i_b/denom) * s * damp
        term = (2.0 * np.pi * i * K / (Z * Z)) * (term_i_b / denom) * s * damp

        c += term

        # relative early-stop
        if np.abs(term) < rel_tol * (1.0 + np.abs(c)):
            break

    return c  # DNT price (PV)


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
    ):
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
        dom_curve: DiscountCurve,
        for_curve: DiscountCurve,
        model,
    ):
        """FX One-Touch Option valuation using the Black-Scholes model
        assuming a continuous (American) barrier from value date to expiry.
        Handles both cash-or-nothing and asset-or-nothing options."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if dom_curve.value_dt != value_dt:
            raise FinError("Domestic Curve date not same as valuation date")

        if for_curve.value_dt != value_dt:
            raise FinError("Foreign Curve date not same as valuation date")

        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        if t_exp <= 0.0:
            raise FinError("Valuation date after or on expiry date.")

        df_d = dom_curve.df(self.expiry_dt)
        df_f = for_curve.df(self.expiry_dt)
        r_d = -np.log(df_d) / t_exp
        r_f = -np.log(df_f) / t_exp

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

        # moved code to Numba for speed
        c = _fast_double_no_touch_pricer(
            S0, L, U, K, t_exp, opt_type_value, r_d, r_f, sigma
        )

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
        dom_curve: DiscountCurve,
        for_curve: DiscountCurve,
        model,
        num_paths: int = 10000,
        num_steps_per_year: int = 252,
        seed: int = 4242,
    ):
        """Double one touch Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Accuracy is not great when compared to the analytical
        result as we only observe the barrier a finite number of times. The
        convergence is slow."""

        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        df_d = dom_curve.df(self.expiry_dt)
        r_d = -np.log(df_d) / t_exp

        df_f = for_curve.df(self.expiry_dt)
        r_f = -np.log(df_f) / t_exp

        num_time_steps = int(t_exp * num_steps_per_year) + 1

        v = model.volatility
        s0 = stock_price
        mu = r_d - r_f

        _, s = get_paths_times(num_paths, num_time_steps, t_exp, mu, s0, v, seed)

        k1 = self.lower_barrier_fx_rate
        k2 = self.upper_barrier_fx_rate

        p_hit = _barrier_pay_at_expiry_double_hit(s, k1, k2)

        if self.option_type == DoubleBarrierTypes.KNOCK_OUT:
            pv = self.payment_size * (1.0 - p_hit) * df_d
        else:
            pv = self.payment_size * p_hit * df_d

        return pv

    ####################################################################################

    def value_mc_bb_slow(
        self,
        value_dt,
        stock_price,
        dom_curve,
        for_curve,
        model,
        num_steps_per_year=52,
        num_paths=100000,
        seed=42,
    ):

        # Double one-touch option valuation with Brownian-bridge correction.
        t_exp = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        if t_exp <= 0.0:
            return 0.0

        df_d = dom_curve.df(self.expiry_dt)
        r_d = -np.log(df_d) / t_exp
        df_f = for_curve.df(self.expiry_dt)
        r_f = -np.log(df_f) / t_exp

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
                    p_hit_lower = np.exp(
                        -2.0 * (x - lnL) * (x_new - lnL) / (sigma**2 * dt)
                    )
                    if np.random.rand() < p_hit_lower:
                        hit = True
                        break

                # Upper barrier
                if x < lnU and x_new < lnU:
                    p_hit_upper = np.exp(
                        -2.0 * (lnU - x) * (lnU - x_new) / (sigma**2 * dt)
                    )
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
        value_dt,
        spot_fx_rate: float,
        dom_curve,
        for_curve,
        model,
        num_steps_per_year: int = 52,
        num_paths: int = 1_000_000,
        seed: int = 42,
    ) -> float:
        """
        PV via Brownian-bridge Monte Carlo with Numba parallelism.
        - KNOCK_IN  (double one-touch): pays K if a barrier is touched before expiry.
        - KNOCK_OUT (double no-touch):   pays K if neither barrier is touched before expiry.
        """
        T = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        if T <= 0.0:
            return 0.0

        df_d = dom_curve.df(self.expiry_dt)
        r_d = -np.log(df_d) / T

        df_f = for_curve.df(self.expiry_dt)
        r_f = -np.log(df_f) / T

        mu = r_d - r_f
        sigma = float(model.volatility)
        S0 = spot_fx_rate
        L = self.lower_barrier_fx_rate
        U = self.upper_barrier_fx_rate

        steps = int(max(1, num_steps_per_year * T))

        # Estimate touch probability with Brownian bridge (parallel)
        p_touch = _p_touch_bb_parallel(S0, L, U, mu, sigma, T, steps, num_paths, seed)

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
