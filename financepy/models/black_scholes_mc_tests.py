##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np

from numba import njit, float64, int64

from math import exp, sqrt

###############################################################################


def value_mc1(s0, t, k, r, q, v, num_paths, seed):
    v_sqrt_t = v * sqrt(t)
    st = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = st * exp(+g[i] * v_sqrt_t)
        payoff += max(s - k, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################


def value_mc2(s0, t, K, r, q, v, num_paths, seed):

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)
    v_sqrt_t = v * np.sqrt(t)
    s = s0 * exp((r - q - v*v / 2.0) * t)

    s = s * np.exp(g * v_sqrt_t)
    payoff = np.maximum(s - K, 0.0)
    average_payoff = np.mean(payoff)

    v = average_payoff * np.exp(-r * t)
    return v

###############################################################################


@njit(float64(float64, float64, float64, float64, float64, float64,
              int64, int64), cache=True, fastmath=True)
def value_mc3(s0, t, K, r, q, v, num_paths, seed):

    v_sqrt_t = v * sqrt(t)
    ss = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = ss * exp(+g[i] * v_sqrt_t)
        payoff += max(s - K, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################


@njit(float64(float64, float64, float64, float64, float64, float64,
              int64, int64), cache=True, fastmath=True)
def value_mc4(s0, t, K, r, q, v, num_paths, seed):

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    v_sqrt_t = v * np.sqrt(t)

    s = s0 * exp((r - q - v*v / 2.0) * t)
    s = s * np.exp(g * v_sqrt_t)

    payoff = np.maximum(s - K, 0.0)
    average_payoff = np.mean(payoff)
    v = average_payoff * np.exp(-r * t)

    return v

###############################################################################


@njit(float64(float64, float64, float64, float64, float64, float64,
              int64, int64), cache=True, fastmath=True)
def value_mc5(s0, t, K, r, q, v, num_paths, seed):

    v_sqrt_t = v * sqrt(t)
    ss = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = ss * exp(+g[i] * v_sqrt_t)
        payoff += max(s - K, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################
