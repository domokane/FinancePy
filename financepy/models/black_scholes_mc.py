# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

########################################################################################
# TODO
########################################################################################

# 1. Verify Sobol generator compatibility with Numba nopython mode.
# 2. Separate random number generation from Monte Carlo pricing kernels.
# 3. Add input validation for t, vol, num_paths and option type.
# 4. Implement control variates using analytic Black-Scholes values.
# 5. Benchmark Sobol, pseudo-random and antithetic convergence.
# 6. Benchmark memory usage of vectorised versus loop-based implementations.
# 7. Improve parallel Monte Carlo reduction structure.
# 8. Add convergence/error analysis utilities.
# 9. Consolidate duplicated Monte Carlo payoff logic.
# 10. Remove unused imports and clean up comments.

########################################################################################

import random
import math
import numpy as np

from numba import njit, prange
from ..utils.global_types import OptionTypes
from ..models.sobol import get_gaussian_sobol
from ..utils.stats import std_err
from ..utils.error import FinError

########################################################################################


@njit(cache=True)
def _validate_mc_inputs(
    t_exp: float,
    k: float,
    option_type_value,
    stock_price: float,
    volatility: float,
    num_path_pairs: int,
) -> None:

    if t_exp < 0.0:
        raise FinError("Time to expiry must be positive.")

    if k < 0.0:
        raise FinError("Strike must be non-negative.")

    if option_type_value != OptionTypes.EUROPEAN_CALL.value and option_type_value != OptionTypes.EUROPEAN_PUT.value:
        raise FinError("Option type must be EUROPEAN call or put.")

    if stock_price <= 0.0:
        raise FinError("Stock price must be positive.")

    if volatility < 0.0:
        raise FinError("Volatility must be non-negative.")

    if num_path_pairs < 2:
        raise FinError("Number of path pairs must be 2 or more.")

########################################################################################


@njit(cache=True, fastmath=True, parallel=False)
def value_at_expiry(
    s: float,
    k: float,
    opt_type: int,
):
    """Value a European option exactly at expiry."""

    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        value = max(s - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        value = max(k - s, 0.0)
    else:
        raise FinError("Option type must be EUROPEAN call or put.")

    return value, 0.0

########################################################################################


def value_mc_nonumba_nonumpy(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    random.seed(seed)

    mu = r - q
    v2 = v * v
    v_sqrt_t = v * math.sqrt(t)

    ss = s * math.exp((mu - 0.5 * v2) * t)

    # Generate Gaussian samples.
    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = None

    sum_payoff = 0.0
    sum_payoff_sq = 0.0

    for i in range(num_paths):

        if use_sobol == 1:
            z = g[i]
        else:
            z = random.gauss(0.0, 1.0)

        m = math.exp(z * v_sqrt_t)

        s_1 = ss * m
        s_2 = ss / m

        if opt_type == OptionTypes.EUROPEAN_CALL.value:
            payoff_1 = max(s_1 - k, 0.0)
            payoff_2 = max(s_2 - k, 0.0)
        else:
            payoff_1 = max(k - s_1, 0.0)
            payoff_2 = max(k - s_2, 0.0)

        payoff = 0.5 * (payoff_1 + payoff_2)

        sum_payoff += payoff
        sum_payoff_sq += payoff * payoff

    mean_payoff = sum_payoff / num_paths

    variance = (sum_payoff_sq - num_paths * mean_payoff * mean_payoff)
    variance = variance / (num_paths - 1)

    # Protect against tiny negative values caused by floating-point rounding.
    variance = max(variance, 0.0)

    error = math.sqrt(variance / num_paths)

    discount = math.exp(-r * t)

    value = discount * mean_payoff
    error = discount * error

    return value, error

########################################################################################


@njit(cache=True, fastmath=True, parallel=False)
def value_mc_numpy_only(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    v_sqrt_t = v * np.sqrt(t)

    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = np.random.standard_normal(num_paths)

    ss = s * np.exp((mu - v2 / 2.0) * t)
    m = np.exp(g * v_sqrt_t)
    s_1 = ss * m
    s_2 = ss / m

    # Not sure if it is correct to do antithetics with sobols but why not ? Well ...
    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoffs_1 = np.maximum(s_1 - k, 0.0)
        payoffs_2 = np.maximum(s_2 - k, 0.0)
    else:
        payoffs_1 = np.maximum(k - s_1, 0.0)
        payoffs_2 = np.maximum(k - s_2, 0.0)

    payoffs = (payoffs_1 + payoffs_2)/2.0

    discount = np.exp(-r * t)
    value = discount * np.mean(payoffs)
    error = discount * std_err(payoffs)

    return value, error

########################################################################################


@njit(cache=True, fastmath=True, parallel=False)
def value_mc_numpy_numba(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    v_sqrt_t = v * np.sqrt(t)

    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = np.random.standard_normal(num_paths)

    ss = s * np.exp((mu - v2 / 2.0) * t)
    m = np.exp(g * v_sqrt_t)
    s_1 = ss * m
    s_2 = ss / m

    # Not sure if it is correct to do antithetics with sobols but why not ?
    if opt_type == OptionTypes.EUROPEAN_CALL.value:
        payoffs_1 = np.maximum(s_1 - k, 0.0)
        payoffs_2 = np.maximum(s_2 - k, 0.0)
    else:
        payoffs_1 = np.maximum(k - s_1, 0.0)
        payoffs_2 = np.maximum(k - s_2, 0.0)

    payoffs = (payoffs_1 + payoffs_2) / 2.0

    discount = np.exp(-r * t)
    value = discount * np.mean(payoffs)
    error = discount * std_err(payoffs)

    return value, error


########################################################################################


@njit(cache=True, fastmath=True, parallel=False)
def value_mc_numba_only(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    random.seed(seed)

    mu = r - q
    v2 = v * v
    v_sqrt_t = v * math.sqrt(t)

    ss = s * math.exp((mu - 0.5 * v2) * t)

    # Generate Gaussian samples.
    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = None

    sum_payoff = 0.0
    sum_payoff_sq = 0.0

    for i in range(num_paths):

        if use_sobol == 1:
            z = g[i]
        else:
            z = random.gauss(0.0, 1.0)

        m = math.exp(z * v_sqrt_t)

        s_1 = ss * m
        s_2 = ss / m

        if opt_type == OptionTypes.EUROPEAN_CALL.value:
            payoff_1 = max(s_1 - k, 0.0)
            payoff_2 = max(s_2 - k, 0.0)
        else:
            payoff_1 = max(k - s_1, 0.0)
            payoff_2 = max(k - s_2, 0.0)

        payoff = 0.5 * (payoff_1 + payoff_2)

        sum_payoff += payoff
        sum_payoff_sq += payoff * payoff

    mean_payoff = sum_payoff / num_paths

    variance = (sum_payoff_sq - num_paths * mean_payoff * mean_payoff)
    variance = variance / (num_paths - 1)

    # Protect against tiny negative values caused by floating-point rounding.
    variance = max(variance, 0.0)

    error = math.sqrt(variance / num_paths)

    discount = math.exp(-r * t)

    value = discount * mean_payoff
    error = discount * error

    return value, error


########################################################################################


@njit(cache=True, fastmath=True, parallel=False)
def value_mc_numba_noanti(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:
    # No use of Numpy vectorisation but NUMBA
    # No use of antithetic variables

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    v_sqrt_t = v * np.sqrt(t)
    payoff = 0.0

    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = np.random.standard_normal(num_paths)

    ss = s * np.exp((mu - v2 / 2.0) * t)

    if opt_type == OptionTypes.EUROPEAN_CALL.value:

        for i in range(0, num_paths):
            gg = g[i]
            s_1 = ss * np.exp(+gg * v_sqrt_t)
            payoff += max(s_1 - k, 0.0)

    else:

        for i in range(0, num_paths):
            gg = g[i]
            s_1 = ss * np.exp(+gg * v_sqrt_t)
            payoff += max(k - s_1, 0.0)

    value = payoff * np.exp(-r * t) / num_paths
    return value, 0.0


########################################################################################


@njit(cache=True, fastmath=True, parallel=True)
def value_mc_numba_parallel(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type: int,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> tuple[float, float]:
    """Black-Scholes Monte Carlo valuation using Numba parallelisation.

    Returns
    -------
    tuple[float, float]
        Option value and Monte Carlo standard error.
    """

    _validate_mc_inputs(
        t,
        k,
        opt_type,
        s,
        v,
        num_paths,
    )

    if t == 0.0:
        return value_at_expiry(s, k, opt_type)

    np.random.seed(seed)

    mu = r - q
    v2 = v**2
    v_sqrt_t = v * np.sqrt(t)

    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = np.random.standard_normal(num_paths)

    ss = s * np.exp((mu - v2 / 2.0) * t)

    # Each entry contains the average payoff of an antithetic pair.
    payoffs = np.empty(num_paths)

    if opt_type == OptionTypes.EUROPEAN_CALL.value:

        for i in prange(num_paths):
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)

            payoff1 = max(s_1 - k, 0.0)
            payoff2 = max(s_2 - k, 0.0)

            payoffs[i] = (payoff1 + payoff2) / 2.0

    else:

        for i in prange(num_paths):
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)

            payoff1 = max(k - s_1, 0.0)
            payoff2 = max(k - s_2, 0.0)

            payoffs[i] = (payoff1 + payoff2) / 2.0

    # Monte Carlo estimate
    average_payoff = np.mean(payoffs)

    discount_factor = np.exp(-r * t)
    value = average_payoff * discount_factor

    # Standard error of the Monte Carlo estimator.
    #
    # Each element of payoffs is one antithetic-pair observation,
    # so the effective sample size here is num_paths.
    payoff_std = np.std(payoffs)
    error = discount_factor * payoff_std / np.sqrt(num_paths)

    return value, error
