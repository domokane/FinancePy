##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np

from numba import njit, float64, int64, prange

from math import exp, sqrt

###############################################################################

def valueMC1(s0, t, K, r, q, v, num_paths, seed):

    vsqrtt = v * sqrt(t)
    ss = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = ss * exp(+g[i] * vsqrtt)
        payoff += max(s - K, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################

def valueMC2(s0, t, K, r, q, v, num_paths, seed):

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)
    vsqrtt = v * np.sqrt(t)
    s = s0 * exp((r - q - v*v / 2.0) * t) 

    s = s * np.exp(g * vsqrtt)
    payoff = np.maximum(s - K, 0.0)
    averagePayoff = np.mean(payoff)

    v = averagePayoff * np.exp(-r * t)
    return v

###############################################################################

@njit(float64(float64, float64, float64, float64, float64, float64, 
             int64, int64), cache=True, fastmath=True)
def valueMC3(s0, t, K, r, q, v, num_paths, seed):

    vsqrtt = v * sqrt(t)
    ss = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in range(0, num_paths):
        s = ss * exp(+g[i] * vsqrtt)
        payoff += max(s - K, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################

@njit(float64(float64, float64, float64, float64, float64, float64, 
             int64, int64), cache=True, fastmath=True)
def valueMC4(s0, t, K, r, q, v, num_paths, seed):

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    vsqrtt = v * np.sqrt(t)

    s = s0 * exp((r - q - v*v / 2.0) * t) 
    s = s * np.exp(g * vsqrtt)

    payoff = np.maximum(s - K, 0.0)
    averagePayoff = np.mean(payoff)
    v = averagePayoff * np.exp(-r * t)

    return v

###############################################################################

@njit(float64(float64, float64, float64, float64, float64, float64, 
             int64, int64), cache=True, fastmath=True, parallel=True)
def valueMC5(s0, t, K, r, q, v, num_paths, seed):

    vsqrtt = v * sqrt(t)
    ss = s0 * exp((r - q - v*v / 2.0) * t)

    np.random.seed(seed)
    g = np.random.standard_normal(num_paths)

    payoff = 0.0
    for i in prange(0, num_paths):
        s = ss * exp(+g[i] * vsqrtt)
        payoff += max(s - K, 0.0)

    v = payoff * np.exp(-r * t) / num_paths

    return v

###############################################################################
