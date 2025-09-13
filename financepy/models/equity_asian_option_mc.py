# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import njit
import numba as nb

from ..utils.error import FinError
from ..utils.math import covar

from ..utils.global_types import OptionTypes

error_str = "In averaging period so need to enter accrued average."


@njit(cache=True, fastmath=True, parallel=False)
def equity_asian_value_mc_numba(
    t0: float,
    t: float,
    tau: float,
    k: float,
    n: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths: int,
    seed: int,
    accrued_average: float,
) -> float:

    if opt_type not in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        raise FinError("Invalid Option Type: Must be EUROPEAN call or put")

    # Start pricing here
    np.random.seed(seed)
    multiplier = 1.0

    if t0 < 0.0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    mu = interest_rate - dividend_yield
    v2 = volatility**2
    dt = (t - t0) / n

    payoff_a = 0.0

    for _ in range(num_paths):

        # evolve stock price to start of averaging period
        g = np.random.standard_normal(1)

        s_1 = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 + g[0] * np.sqrt(t0) * volatility
        )
        s_2 = stock_price * np.exp(
            (mu - v2 / 2.0) * t0 - g[0] * np.sqrt(t0) * volatility
        )

        # enter averaging period
        s_1_arithmetic = 0.0
        s_2_arithmetic = 0.0

        g = np.random.standard_normal(n)

        exp_drift = np.exp((mu - v2 / 2.0) * dt)

        for obs in range(0, n):

            s_1 = s_1 * exp_drift * np.exp(+g[obs] * np.sqrt(dt) * volatility)
            s_2 = s_2 * exp_drift * np.exp(-g[obs] * np.sqrt(dt) * volatility)

            s_1_arithmetic += s_1
            s_2_arithmetic += s_2

        s_1_arithmetic /= n
        s_2_arithmetic /= n

        if opt_type == OptionTypes.EUROPEAN_CALL:
            payoff_a += max(s_1_arithmetic - k, 0.0)
            payoff_a += max(s_2_arithmetic - k, 0.0)
        else:
            payoff_a += max(k - s_1_arithmetic, 0.0)
            payoff_a += max(k - s_2_arithmetic, 0.0)

    v_a = payoff_a * np.exp(-interest_rate * t) / num_paths / 2.0
    v_a = v_a * multiplier
    return v_a


########################################################################################


@njit(cache=True, parallel=False)
def equity_asian_value_mc_fast_numba(
    t0: float,
    t: float,
    tau: float,
    k: float,
    n: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths: int,
    seed: int,
    accrued_average: float,
) -> float:

    np.random.seed(seed)
    mu = interest_rate - dividend_yield
    s0 = stock_price
    v2 = volatility**2
    dt = (t - t0) / n
    r = interest_rate
    num_paths = int(num_paths)

    multiplier = 1.0

    if t0 < 0.0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    # g = np.random.normal(0.0, 1.0, size=(num_paths))

    gg = np.empty(num_paths, np.float64)

    #    for ip in nb.prange(num_paths):
    for ip in range(0, num_paths):
        rv = np.random.normal()
        gg[ip] = rv

    s_1 = np.empty(num_paths)
    s_2 = np.empty(num_paths)

    exp_drift_t0 = np.exp((mu - v2 / 2.0) * t0)

    for ip in range(0, num_paths):
        #    for ip in nb.prange(num_paths):
        s_1[ip] = s0 * exp_drift_t0 * np.exp(+gg[ip] * np.sqrt(t0) * volatility)
        s_2[ip] = s0 * exp_drift_t0 * np.exp(-gg[ip] * np.sqrt(t0) * volatility)

    s_1_arithmetic = np.zeros(num_paths)
    s_2_arithmetic = np.zeros(num_paths)
    sigma_root_dt = volatility * np.sqrt(dt)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths)

        #        for ip in nb.prange(num_paths):
        for ip in range(0, num_paths):
            s_1[ip] = s_1[ip] * np.exp((mu - v2 / 2.0) * dt + g[ip] * sigma_root_dt)
            s_2[ip] = s_2[ip] * np.exp((mu - v2 / 2.0) * dt - g[ip] * sigma_root_dt)

        for ip in nb.prange(num_paths):
            s_1_arithmetic[ip] += s_1[ip] / n
            s_2_arithmetic[ip] += s_2[ip] / n

    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
    else:
        raise FinError("Unknown option type.")

    payoff_a = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    v_a = multiplier * payoff_a * np.exp(-r * t) / 2.0
    return v_a


########################################################################################


@njit(cache=True, parallel=False)
def equity_asian_value_mc_fast_cv_numba(
    t0: float,
    t: float,
    tau: float,
    k: float,
    n: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths: int,
    seed: int,
    accrued_average: float,
    v_g_exact: float,
) -> float:

    np.random.seed(seed)
    mu = interest_rate - dividend_yield
    v2 = volatility**2
    dt = (t - t0) / n
    r = interest_rate

    multiplier = 1.0

    if t0 < 0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # we adjust the strike to account for the accrued coupon
        k = (k * tau + accrued_average * t0) / t
        # the number of options is rescaled also
        multiplier = t / tau
        # there is no pre-averaging time
        t0 = 0.0
        # the number of observations is scaled and floored at 1
        n = int(n * t / tau + 0.5) + 1

    # evolve stock price to start of averaging period
    g = np.random.normal(0.0, 1.0, size=num_paths)

    s_1 = np.empty(num_paths)
    s_2 = np.empty(num_paths)

    sigma_root_t0 = np.sqrt(t0) * volatility

    for ip in range(0, num_paths):
        #    for ip in nb.prange(num_paths):

        s_1[ip] = stock_price * np.exp((mu - v2 / 2.0) * t0 + g[ip] * sigma_root_t0)
        s_2[ip] = stock_price * np.exp((mu - v2 / 2.0) * t0 - g[ip] * sigma_root_t0)

    s_1_arithmetic = np.zeros(num_paths)
    s_2_arithmetic = np.zeros(num_paths)
    ln_s_1_geometric = np.zeros(num_paths)
    ln_s_2_geometric = np.zeros(num_paths)

    sigma_root_dt = np.sqrt(dt) * volatility

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths)

        for ip in range(0, num_paths):
            #        for ip in nb.prange(num_paths):

            s_1[ip] = s_1[ip] * np.exp((mu - v2 / 2.0) * dt + g[ip] * sigma_root_dt)
            s_2[ip] = s_2[ip] * np.exp((mu - v2 / 2.0) * dt - g[ip] * sigma_root_dt)

        for ip in nb.prange(num_paths):
            s_1_arithmetic[ip] += s_1[ip]
            s_2_arithmetic[ip] += s_2[ip]
            ln_s_1_geometric[ip] += np.log(s_1[ip])
            ln_s_2_geometric[ip] += np.log(s_2[ip])

    s_1_geometric = np.empty(num_paths)
    s_2_geometric = np.empty(num_paths)

    for ip in nb.prange(num_paths):
        s_1_arithmetic[ip] /= n
        s_1_geometric[ip] = np.exp(ln_s_1_geometric[ip] / n)
        s_2_arithmetic[ip] /= n
        s_2_geometric[ip] = np.exp(ln_s_2_geometric[ip] / n)

    if opt_type == OptionTypes.EUROPEAN_CALL:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_g_1 = np.maximum(s_1_geometric - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
        payoff_g_2 = np.maximum(s_2_geometric - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_g_1 = np.maximum(k - s_1_geometric, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
        payoff_g_2 = np.maximum(k - s_2_geometric, 0.0)
    else:
        raise FinError("Unknown Option Type")

    payoff_a = np.concatenate((payoff_a_1, payoff_a_2), axis=0)
    payoff_g = np.concatenate((payoff_g_1, payoff_g_2), axis=0)

    # Now we do the control variate adjustment
    m = covar(payoff_a, payoff_a)

    if np.abs(m[1][1]) < 1e-10:
        lam = 0.0
    else:
        lam = m[0][1] / m[1][1]

    payoff_a_mean = np.mean(payoff_a)
    payoff_g_mean = np.mean(payoff_g)

    v_a = payoff_a_mean * np.exp(-r * t) * multiplier
    v_g = payoff_g_mean * np.exp(-r * t) * multiplier

    epsilon = v_g_exact - v_g
    v_a_cv = v_a + lam * epsilon

    return v_a_cv
