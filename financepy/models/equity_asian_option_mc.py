# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import njit

from ..utils.error import FinError
from ..utils.math import covar

from ..utils.global_types import OptionTypes

error_str = "In averaging period so need to enter accrued average."


@njit(cache=True)
def _validate_asian_mc_inputs(
    t_avg: float,
    t_exp: float,
    k: float,
    num_obs_per_year: int,
    stock_price: float,
    volatility: float,
    num_paths: int,
) -> None:

    if t_exp <= 0.0:
        raise FinError("Time to expiry must be positive.")

    if t_avg >= t_exp:
        raise FinError("Averaging start time t0 must be less than expiry t.")

    if k < 0.0:
        raise FinError("Strike must be non-negative.")

    if num_obs_per_year <= 0:
        raise FinError("Number of observations per year must be positive.")

    if stock_price <= 0.0:
        raise FinError("Stock price must be positive.")

    if volatility < 0.0:
        raise FinError("Volatility must be non-negative.")

    if num_paths <= 0:
        raise FinError("Number of paths must be positive.")


################################################################################

@njit(cache=True, fastmath=True, parallel=False)
def equity_asian_value_mc_numba(
    t_avg: float,
    t_exp: float,
    k: float,
    num_obs_per_year: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths_over_two: int,
    seed: int,
    accrued_average: float,
) -> float:

    if opt_type not in [
        OptionTypes.EUROPEAN_CALL.value,
        OptionTypes.EUROPEAN_PUT.value,
    ]:
        raise FinError("Invalid Option Type: Must be EUROPEAN call or put")

    tau = t_exp - t_avg

    _validate_asian_mc_inputs(
        t_avg,
        t_exp,
        k,
        num_obs_per_year,
        stock_price,
        volatility,
        num_paths_over_two,
    )

    # Start pricing here
    np.random.seed(seed)
    multiplier = 1.0

    if t_avg < 0.0:

        if accrued_average is None:
            raise FinError(error_str)

        k = (k * tau + accrued_average * t_avg) / t_exp
        multiplier = t_exp / tau
        t_avg = 0.0

    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    dt = averaging_time / n

    mu = interest_rate - dividend_yield
    v2 = volatility**2

    payoff_a = 0.0

    for _ in range(num_paths_over_two):

        # evolve stock price to start of averaging period
        g = np.random.normal()

        s_1 = stock_price * np.exp((mu - v2 / 2.0) * t_avg + g * np.sqrt(t_avg) * volatility)
        s_2 = stock_price * np.exp((mu - v2 / 2.0) * t_avg - g * np.sqrt(t_avg) * volatility)

        # enter averaging period
        s_1_arithmetic = 0.0
        s_2_arithmetic = 0.0

        g = np.random.normal(0.0, 1.0, size=n)

        exp_drift = np.exp((mu - v2 / 2.0) * dt)

        for obs in range(0, n):

            s_1 = s_1 * exp_drift * np.exp(+g[obs] * np.sqrt(dt) * volatility)
            s_2 = s_2 * exp_drift * np.exp(-g[obs] * np.sqrt(dt) * volatility)

            s_1_arithmetic += s_1
            s_2_arithmetic += s_2

        s_1_arithmetic /= n
        s_2_arithmetic /= n

        if opt_type == OptionTypes.EUROPEAN_CALL.value:
            payoff_a += max(s_1_arithmetic - k, 0.0)
            payoff_a += max(s_2_arithmetic - k, 0.0)
        elif opt_type == OptionTypes.EUROPEAN_PUT.value:
            payoff_a += max(k - s_1_arithmetic, 0.0)
            payoff_a += max(k - s_2_arithmetic, 0.0)
        else:
            return np.nan

    v_a = payoff_a * np.exp(-interest_rate * t_exp) / num_paths_over_two / 2.0
    v_a = v_a * multiplier
    return v_a


########################################################################################


@njit(cache=True, parallel=False)
def equity_asian_value_mc_fast_numba(
    t_avg: float,
    t_exp: float,
    k: float,
    num_obs_per_year: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths_over_two: int,
    seed: int,
    accrued_average: float,
) -> float:

    _validate_asian_mc_inputs(
        t_avg,
        t_exp,
        k,
        num_obs_per_year,
        stock_price,
        volatility,
        num_paths_over_two,
    )

    np.random.seed(seed)
    mu = interest_rate - dividend_yield
    s0 = stock_price

    v2 = volatility**2
    r = interest_rate
    num_paths_over_two = int(num_paths_over_two)

    multiplier = 1.0
    tau = t_exp - t_avg

    if t_avg < 0.0:

        if accrued_average is None:
            raise FinError(error_str)

        # Adjust strike for observations already accrued.
        k = (k * tau + accrued_average * t_avg) / t_exp

        # Rescale the option value.
        multiplier = t_exp / tau

        # Remaining averaging period starts now.
        t_avg = 0.0

    # Number of future observations.
    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    # Ensure the final observation is at expiry.
    dt = averaging_time / n

    gg = np.empty(num_paths_over_two, np.float64)

    for ip in range(0, num_paths_over_two):
        rv = np.random.normal()
        gg[ip] = rv

    s_1 = np.empty(num_paths_over_two)
    s_2 = np.empty(num_paths_over_two)

    exp_drift_t_avg = np.exp((mu - v2 / 2.0) * t_avg)

    for ip in range(0, num_paths_over_two):
        s_1[ip] = s0 * exp_drift_t_avg * np.exp(+gg[ip] * np.sqrt(t_avg) * volatility)
        s_2[ip] = s0 * exp_drift_t_avg * np.exp(-gg[ip] * np.sqrt(t_avg) * volatility)

    s_1_arithmetic = np.zeros(num_paths_over_two)
    s_2_arithmetic = np.zeros(num_paths_over_two)
    sigma_root_dt = volatility * np.sqrt(dt)

    exp_drift_dt = np.exp((mu - v2 / 2.0) * dt)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths_over_two)

        for ip in range(0, num_paths_over_two):
            s_1[ip] = s_1[ip] * exp_drift_dt * np.exp(+g[ip] * sigma_root_dt)
            s_2[ip] = s_2[ip] * exp_drift_dt * np.exp(-g[ip] * sigma_root_dt)

        for ip in range(0, num_paths_over_two):
            s_1_arithmetic[ip] += s_1[ip] / n
            s_2_arithmetic[ip] += s_2[ip] / n

    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
    else:
        return np.nan

    payoff_a = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    v_a = multiplier * payoff_a * np.exp(-r * t_exp) / 2.0
    return v_a


########################################################################################


@njit(cache=True, parallel=False)
def equity_asian_value_mc_fast_cv_numba(
    t_avg: float,
    t_exp: float,
    k: float,
    num_obs_per_year: int,
    opt_type: int,
    stock_price: float,
    interest_rate: float,
    dividend_yield: float,
    volatility: float,
    num_paths_over_two: int,
    seed: int,
    accrued_average: float,
    v_g_exact: float,
) -> float:

    _validate_asian_mc_inputs(
        t_avg,
        t_exp,
        k,
        num_obs_per_year,
        stock_price,
        volatility,
        num_paths_over_two,
    )

    np.random.seed(seed)

    mu = interest_rate - dividend_yield
    v2 = volatility**2
    r = interest_rate
    tau = t_exp - t_avg

    multiplier = 1.0

    if t_avg < 0.0:  # we are in the averaging period

        if accrued_average is None:
            raise FinError(error_str)

        # Adjust strike for the accrued part of the average.
        k = (k * tau + accrued_average * t_avg) / t_exp

        # Rescale the number of options.
        multiplier = t_exp / tau

        # Remaining averaging period starts now.
        t_avg = 0.0

    # Number of future observations over the remaining averaging period.
    averaging_time = t_exp - t_avg

    n = max(
        1,
        int(averaging_time * num_obs_per_year + 0.5),
    )

    # Make the final observation fall exactly at expiry.
    dt = averaging_time / n

    # evolve stock price to start of averaging period
    g = np.random.normal(0.0, 1.0, size=num_paths_over_two)

    s_1 = np.empty(num_paths_over_two)
    s_2 = np.empty(num_paths_over_two)

    sigma_root_t0 = np.sqrt(t_avg) * volatility
    s_exp_drift_t0 = stock_price * np.exp((mu - v2 / 2.0) * t_avg)

    for ip in range(0, num_paths_over_two):
        s_1[ip] = s_exp_drift_t0 * np.exp(+g[ip] * sigma_root_t0)
        s_2[ip] = s_exp_drift_t0 * np.exp(-g[ip] * sigma_root_t0)

    s_1_arithmetic = np.zeros(num_paths_over_two)
    s_2_arithmetic = np.zeros(num_paths_over_two)
    ln_s_1_geometric = np.zeros(num_paths_over_two)
    ln_s_2_geometric = np.zeros(num_paths_over_two)

    sigma_root_dt = np.sqrt(dt) * volatility
    exp_drift_dt = np.exp((mu - v2 / 2.0) * dt)

    for _ in range(0, n):

        g = np.random.normal(0.0, 1.0, size=num_paths_over_two)

        for ip in range(0, num_paths_over_two):
            s_1[ip] = s_1[ip] * exp_drift_dt * np.exp(+g[ip] * sigma_root_dt)
            s_2[ip] = s_2[ip] * exp_drift_dt * np.exp(-g[ip] * sigma_root_dt)

        for ip in range(num_paths_over_two):
            s_1_arithmetic[ip] += s_1[ip]
            s_2_arithmetic[ip] += s_2[ip]
            ln_s_1_geometric[ip] += np.log(s_1[ip])
            ln_s_2_geometric[ip] += np.log(s_2[ip])

    s_1_geometric = np.empty(num_paths_over_two)
    s_2_geometric = np.empty(num_paths_over_two)

    for ip in range(num_paths_over_two):
        s_1_arithmetic[ip] /= n
        s_1_geometric[ip] = np.exp(ln_s_1_geometric[ip] / n)
        s_2_arithmetic[ip] /= n
        s_2_geometric[ip] = np.exp(ln_s_2_geometric[ip] / n)

    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoff_a_1 = np.maximum(s_1_arithmetic - k, 0.0)
        payoff_g_1 = np.maximum(s_1_geometric - k, 0.0)
        payoff_a_2 = np.maximum(s_2_arithmetic - k, 0.0)
        payoff_g_2 = np.maximum(s_2_geometric - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1_arithmetic, 0.0)
        payoff_g_1 = np.maximum(k - s_1_geometric, 0.0)
        payoff_a_2 = np.maximum(k - s_2_arithmetic, 0.0)
        payoff_g_2 = np.maximum(k - s_2_geometric, 0.0)
    else:
        return np.nan

    payoff_a = np.concatenate((payoff_a_1, payoff_a_2), axis=0)
    payoff_g = np.concatenate((payoff_g_1, payoff_g_2), axis=0)

    # Now we do the control variate adjustment
    m = covar(payoff_a, payoff_g)

    if np.abs(m[1][1]) < 1e-10:
        lam = 0.0
    else:
        lam = m[0][1] / m[1][1]

    payoff_a_mean = np.mean(payoff_a)
    payoff_g_mean = np.mean(payoff_g)

    v_a = payoff_a_mean * np.exp(-r * t_exp) * multiplier
    v_g = payoff_g_mean * np.exp(-r * t_exp) * multiplier

    epsilon = v_g_exact - v_g
    v_a_cv = v_a + lam * epsilon

    return v_a_cv
