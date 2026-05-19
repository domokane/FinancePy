# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Use for Numba timing comparisons

import numpy as np
from numba import njit, float64, int64

########################################################################################


def value_mc1(
    s0: float,
    t: float,
    k: float,
    r: float,
    q: float,
    vol: float,
    num_paths: int,
    seed: int,
) -> float:

    v_sqrt_t = vol * np.sqrt(t)
    st = s0 * np.exp((r - q - vol * vol / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = st * np.exp(+g[i] * v_sqrt_t)
        payoff += np.maximum(s - k, 0.0)

    value = payoff * np.exp(-r * t) / num_paths
    return value


########################################################################################


def value_mc2(
    s0: float,
    t: float,
    k: float,
    r: float,
    q: float,
    vol: float,
    num_paths: int,
    seed: int,
) -> float:

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)
    v_sqrt_t = vol * np.sqrt(t)
    s = s0 * np.exp((r - q - vol * vol / 2.0) * t)
    s = s * np.exp(g * v_sqrt_t)
    payoff = np.maximum(s - k, 0.0)
    average_payoff = np.mean(payoff)
    value = average_payoff * np.exp(-r * t)
    return value


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64, float64, int64, int64),
    cache=True,
    fastmath=True,
)
def value_mc3(
    s0: float,
    t: float,
    k: float,
    r: float,
    q: float,
    vol: float,
    num_paths: int,
    seed: int,
) -> float:

    v_sqrt_t = vol * np.sqrt(t)
    ss = s0 * np.exp((r - q - vol * vol / 2.0) * t)
    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = ss * np.exp(+g[i] * v_sqrt_t)
        payoff += np.maximum(s - k, 0.0)

    value = payoff * np.exp(-r * t) / num_paths
    return value


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64, float64, int64, int64),
    cache=True,
    fastmath=True,
)
def value_mc4(
    s0: float,
    t: float,
    k: float,
    r: float,
    q: float,
    vol: float,
    num_paths: int,
    seed: int,
) -> float:

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)
    v_sqrt_t = vol * np.sqrt(t)
    s = s0 * np.exp((r - q - vol * vol / 2.0) * t)
    s = s * np.exp(g * v_sqrt_t)
    payoff = np.maximum(s - k, 0.0)
    average_payoff = np.mean(payoff)
    value = average_payoff * np.exp(-r * t)
    return value


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64, float64, int64, int64),
    cache=True,
    fastmath=True,
)
def value_mc5(s0, t, k, r, q, vol, num_paths, seed):

    np.random.seed(seed)

    if num_paths % 2 == 1:
        num_paths -= 1

    half_paths = num_paths // 2
    g = np.random.standard_normal(half_paths)
    g = np.concatenate((g, -g))

    v_sqrt_t = vol * np.sqrt(t)
    s = s0 * np.exp((r - q - vol * vol / 2.0) * t)
    s = s * np.exp(g * v_sqrt_t)

    payoff = np.maximum(s - k, 0.0)
    value = np.mean(payoff) * np.exp(-r * t)

    return value


########################################################################################
