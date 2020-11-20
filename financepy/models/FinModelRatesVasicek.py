##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, exp
from numba import njit, float64, int64
import numpy as np

from ..finutils.FinHelperFunctions import labelToString

##########################################################################
# dr = a(b-r) + sigma dW
##########################################################################

# TO DO - DECIDE WHETHER TO OO MODEL

###############################################################################


class FinModelRatesVasicek():

    def __init__(self, a, b, sigma):
        self._a = a
        self._b = b
        self._sigma = sigma

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("a", self._a)
        s += labelToString("b", self._b)
        s += labelToString("sigma", self._sigma)
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
def zeroPrice(r0, a, b, sigma, t):
    B = (1.0 - exp(-a * t)) / a
    A = exp((b - sigma * sigma / 2.0 / a / a) *
            (B - t) - B * B * sigma * sigma / 4.0 / a)
    zcb = A * exp(-r0 * B)
    return zcb

###############################################################################


@njit(float64[:](float64, float64, float64, float64, float64, float64, int64))
def ratePath_MC(r0, a, b, sigma, t, dt, seed):

    np.random.seed(seed)
    numSteps = int(t / dt)
    ratePath = np.zeros(numSteps)
    ratePath[0] = r0
    numPaths = 1

    sigmasqrtdt = sigma * sqrt(dt)

    for iPath in range(0, numPaths):

        r = r0
        z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

        for iStep in range(1, numSteps):
            r = r + a * (b - r) * dt + z[iStep - 1] * sigmasqrtdt
            ratePath[iStep] = r

    return ratePath

###############################################################################


@njit(float64(float64, float64, float64, float64, float64,
      float64, int64, int64), fastmath=True, cache=True)
def zeroPrice_MC(r0, a, b, sigma, t, dt, numPaths, seed):

    np.random.seed(seed)
    numSteps = int(t / dt)
    sigmasqrtdt = sigma * sqrt(dt)
    zcb = 0.0
    for iPath in range(0, numPaths):
        z = np.random.normal(0.0, 1.0, size=(numSteps))
        rsum = 0.0
        r = r0
        for iStep in range(0, numSteps):
            r += a * (b - r) * dt + z[iStep] * sigmasqrtdt
            rsum += r * dt
        zcb += exp(-rsum)
    zcb /= numPaths
    return zcb

###############################################################################
