##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, exp
from numba import njit, float64, int64
import numpy as np

from ..utils.helpers import label_to_string

##########################################################################
# dr = a(b-r) + sigma dW
##########################################################################

# TO DO - DECIDE WHETHER TO OO MODEL

###############################################################################


class ModelRatesVasicek():

    def __init__(self, a, b, sigma):
        self._a = a
        self._b = b
        self._sigma = sigma

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("a", self._a)
        s += label_to_string("b", self._b)
        s += label_to_string("sigma", self._sigma)
        return s

###############################################################################


@njit(fastmath=True, cache=True)
def meanr(r0, a, b, t):
    mr = r0 * exp(-a * t) + b * (1 - exp(-a * t))
    return mr

###############################################################################


@njit(fastmath=True, cache=True)
def variancer(a, b, sigma, t):
    vr = sigma * sigma * (1.0 - exp(-2.0 * a * t)) / 2.0 / a
    return vr

###############################################################################


@njit(fastmath=True, cache=True)
def zero_price(r0, a, b, sigma, t):
    B = (1.0 - exp(-a * t)) / a
    A = exp((b - sigma * sigma / 2.0 / a / a) *
            (B - t) - B * B * sigma * sigma / 4.0 / a)
    zcb = A * exp(-r0 * B)
    return zcb

###############################################################################


@njit(float64[:](float64, float64, float64, float64, float64, float64, int64))
def rate_path_mc(r0, a, b, sigma, t, dt, seed):

    np.random.seed(seed)
    num_steps = int(t / dt)
    rate_path = np.zeros(num_steps)
    rate_path[0] = r0
    num_paths = 1

    sigmasqrt_dt = sigma * sqrt(dt)

    for iPath in range(0, num_paths):

        r = r0
        z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

        for iStep in range(1, num_steps):
            r = r + a * (b - r) * dt + z[iStep - 1] * sigmasqrt_dt
            rate_path[iStep] = r

    return rate_path

###############################################################################


@njit(float64(float64, float64, float64, float64, float64,
              float64, int64, int64), fastmath=True, cache=True)
def zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed):

    np.random.seed(seed)
    num_steps = int(t / dt)
    sigmasqrt_dt = sigma * sqrt(dt)
    zcb = 0.0
    for iPath in range(0, num_paths):
        z = np.random.normal(0.0, 1.0, size=(num_steps))
        rsum = 0.0
        r = r0
        for iStep in range(0, num_steps):
            r += a * (b - r) * dt + z[iStep] * sigmasqrt_dt
            rsum += r * dt
        zcb += exp(-rsum)
    zcb /= num_paths
    return zcb

###############################################################################
