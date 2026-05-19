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

import numpy as np

from numba import njit, float64, int64, prange
from ..utils.global_types import OptionTypes
from ..models.sobol import get_gaussian_sobol

########################################################################################


def _value_mc_nonumba_nonumpy(
    s: float,
    t: float,
    k: float,
    opt_type: int,
    r: float,
    q: float,
    v: float,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> float:
    # SLOWEST - No use of NUMPY vectorisation or NUMBA

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
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)
            payoff += max(s_1 - k, 0.0)
            payoff += max(s_2 - k, 0.0)

    elif opt_type == OptionTypes.EUROPEAN_PUT.value:

        for i in range(0, num_paths):
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)
            payoff += max(k - s_1, 0.0)
            payoff += max(k - s_2, 0.0)

    else:
        return np.nan

    value = payoff * np.exp(-r * t) / num_paths / 2.0
    return value


########################################################################################


def _value_mc_numpy_only(
    s: float,
    t: float,
    k: float,
    opt_type: int,
    r: float,
    q: float,
    v: float,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> float:
    # Use of NUMPY ONLY

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
        payoff_a_1 = np.maximum(s_1 - k, 0.0)
        payoff_a_2 = np.maximum(s_2 - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1, 0.0)
        payoff_a_2 = np.maximum(k - s_2, 0.0)
    else:
        return np.nan

    payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    value = payoff * np.exp(-r * t) / 2.0
    return value


########################################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        int64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64,
    ),
    cache=True,
    fastmath=True,
)
def _value_mc_numpy_numba(
    s: float,
    t: float,
    k: float,
    opt_type: int,
    r: float,
    q: float,
    v: float,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> float:
    # Use of NUMPY ONLY

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
        payoff_a_1 = np.maximum(s_1 - k, 0.0)
        payoff_a_2 = np.maximum(s_2 - k, 0.0)
    elif opt_type == OptionTypes.EUROPEAN_PUT.value:
        payoff_a_1 = np.maximum(k - s_1, 0.0)
        payoff_a_2 = np.maximum(k - s_2, 0.0)
    else:
        return np.nan

    payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
    value = payoff * np.exp(-r * t) / 2.0
    return value


########################################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        int64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64,
    ),
    fastmath=True,
    cache=True,
)
def _value_mc_numba_only(
    s: float,
    t: float,
    k: float,
    opt_type: int,
    r: float,
    q: float,
    v: float,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> float:
    # No use of Numpy vectorisation but NUMBA

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
            s_2 = ss * np.exp(-gg * v_sqrt_t)
            payoff += max(s_1 - k, 0.0)
            payoff += max(s_2 - k, 0.0)

    elif opt_type == OptionTypes.EUROPEAN_PUT.value:

        for i in range(0, num_paths):
            gg = g[i]
            s_1 = ss * np.exp(+gg * v_sqrt_t)
            s_2 = ss * np.exp(-gg * v_sqrt_t)
            payoff += max(k - s_1, 0.0)
            payoff += max(k - s_2, 0.0)

    else:
        return np.nan

    value = payoff * np.exp(-r * t) / num_paths / 2.0
    return value


########################################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        int64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64,
    ),
    fastmath=True,
    cache=True,
    parallel=True,
)
def _value_mc_numba_parallel(
    s: float,
    t: float,
    k: float,
    opt_type: int,
    r: float,
    q: float,
    v: float,
    num_paths: int,
    seed: int,
    use_sobol: int,
) -> float:
    # No use of Numpy vectorisation but NUMBA

    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    v_sqrt_t = v * np.sqrt(t)

    if use_sobol == 1:
        g = get_gaussian_sobol(num_paths, 1)[:, 0]
    else:
        g = np.random.standard_normal(num_paths)

    ss = s * np.exp((mu - v2 / 2.0) * t)

    payoffs = np.empty(num_paths)

    if opt_type == OptionTypes.EUROPEAN_CALL.value:

        for i in prange(0, num_paths):
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)
            payoff1 = max(s_1 - k, 0.0)
            payoff2 = max(s_2 - k, 0.0)
            path_payoff = (payoff1 + payoff2) / 2.0
            payoffs[i] = path_payoff

    elif opt_type == OptionTypes.EUROPEAN_PUT.value:

        for i in prange(0, num_paths):
            s_1 = ss * np.exp(+g[i] * v_sqrt_t)
            s_2 = ss * np.exp(-g[i] * v_sqrt_t)
            payoff1 = max(k - s_1, 0.0)
            payoff2 = max(k - s_2, 0.0)
            path_payoff = (payoff1 + payoff2) / 2.0
            payoffs[i] = path_payoff

    else:
        return np.nan

    average_payoff = np.mean(payoffs)
    value = average_payoff * np.exp(-r * t)
    return value
