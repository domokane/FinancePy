##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np
from ..utils.helpers import label_to_string

###############################################################################
# CIR Process
# dr = a(b-r) + sigma sqrt(r) dW
# Note that r can hit zero if 2.0 * a * b < sigma*sigma:
###############################################################################

# TO DO - DECIDE WHETHER TO OO MODEL
# CAN DO Z SCALING INSIDE NUMPY ?
# ANTITHETICS

from enum import Enum


class CIRNumericalScheme(Enum):
    EULER = 1
    LOGNORMAL = 2
    MILSTEIN = 3
    KAHLJACKEL = 4
    EXACT = 5  # SAMPLES EXACT DISTRIBUTION


###############################################################################

# THIS CLASS IS NOT USED BUT MAY BE USED IF WE CREATE AN OO FRAMEWORK


class CIR_MC():
    ''' This is a Monte Carlo implementation of the Cox-Ingersoll-Ross Model'''

    def __init__(self, a, b, sigma):
        self._a = a
        self._b = b
        self._sigma = sigma

    def __repr__(self):
        """ Return string with class details. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Sigma", self._sigma)
        s += label_to_string("a", self._a)
        s += label_to_string("b", self._b)
        return s


###############################################################################

@njit(fastmath=True, cache=True)
def meanr(r0, a, b, t):
    """ Mean value of a CIR process after time t """
    mr = r0 * np.exp(-a * t) + b * (1.0 - np.exp(-a * t))
    return mr


###############################################################################


@njit(fastmath=True, cache=True)
def variancer(r0, a, b, sigma, t):
    """ Variance of a CIR process after time t """
    vr = r0 * sigma * sigma * (np.exp(-a * t) - np.exp(-2.0 * a * t)) / a
    vr += b * sigma * sigma * ((1.0 - np.exp(-a * t)) ** 2) / 2.0 / a
    return vr


###############################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64),
    fastmath=True,
    cache=True)
def zero_price(r0, a, b, sigma, t):
    """ Price of a zero coupon bond in CIR model. """
    h = np.sqrt(a * a + 2.0 * sigma * sigma)
    denom = 2.0 * h + (a + h) * (np.exp(h * t) - 1.0)
    A = (2.0 * h * np.exp((a + h) * t / 2.0) /
         denom) ** (2.0 * a * b / sigma / sigma)
    B = 2.0 * (np.exp(h * t) - 1.0) / denom
    zcb = A * np.exp(-r0 * B)
    return zcb


###############################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64),
    fastmath=True,
    cache=True)
def draw(rt, a, b, sigma, dt):
    """ Draw a next rate from the CIR model in Monte Carlo. """
    sigma2 = sigma * sigma
    d = 4.0 * a * b / sigma2
    ll = 4.0 * a * np.exp(-a * dt) / sigma2 / (1.0 - np.exp(-a * dt)) * rt
    c = sigma2 * (1.0 - np.exp(-a * dt)) / 4.0 / a

    if d > 1:
        Z = np.random.normal()
        X = np.random.chisquare(d - 1)
        r = c * (X + (Z + np.sqrt(ll)) ** 2)
    else:
        N = np.random.poisson(ll / 2.0)
        X = np.random.chisquare(d + 2 * N)
        r = c * X

    return r


###############################################################################


@njit(
    float64[:](
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        int64,
        int64))
def rate_path_mc(r0, a, b, sigma, t, dt, seed, scheme):
    """ Generate a path of CIR rates using a number of numerical schemes. """

    np.random.seed(seed)
    num_steps = int(t / dt)
    rate_path = np.zeros(num_steps)
    rate_path[0] = r0
    num_paths = 1

    if scheme == CIRNumericalScheme.EULER.value:

        sigmasqrt_dt = sigma * np.sqrt(dt)

        for iPath in range(0, num_paths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrt_dt * np.sqrt(max(r, 0.0))
                rate_path[iStep] = r

    elif scheme == CIRNumericalScheme.LOGNORMAL.value:

        x = np.exp(-a * dt)
        y = 1.0 - x

        for iPath in range(0, num_paths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                mean = x * r + b * y
                var = sigma * sigma * y * (x * r + 0.50 * b * y) / a
                sig = np.sqrt(np.log(1.0 + var / (mean * mean)))
                r = mean * np.exp(-0.5 * sig * sig + sig * z[iStep - 1])
                rate_path[iStep] = r

    elif scheme == CIRNumericalScheme.MILSTEIN.value:

        sigmasqrt_dt = sigma * np.sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0

        for iPath in range(0, num_paths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrt_dt * np.sqrt(max(0.0, r))
                r = r + sigma2dt * (z[iStep - 1] ** 2 - 1.0)
                rate_path[iStep] = r

    elif scheme == CIRNumericalScheme.KAHLJACKEL.value:

        bhat = b - sigma * sigma / 4.0 / a
        sqrt_dt = np.sqrt(dt)

        for iPath in range(0, num_paths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                beta = z[iStep - 1] / sqrt_dt
                rootr = np.sqrt(max(r, 1e-8))
                c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt / 4.0 / rootr
                r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt
                rate_path[iStep] = r

    elif scheme == CIRNumericalScheme.EXACT.value:

        for iPath in range(0, num_paths):

            r = r0

            for iStep in range(1, num_steps):
                r = draw(r, a, b, sigma, dt)
                rate_path[iStep] = r

    return rate_path


###############################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64))
def zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed, scheme):
    """ Determine the CIR zero price using Monte Carlo. """

    if t == 0.0:
        return 1.0

    np.random.seed(seed)

    num_steps = int(t / dt)
    zcb = 0.0

    if scheme == CIRNumericalScheme.EULER.value:

        sigmasqrt_dt = sigma * np.sqrt(dt)

        for iPath in range(0, num_paths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                r_prev = r
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrt_dt * np.sqrt(max(r, 0.0))
                rsum += (r + r_prev)

            zcb += np.exp(-0.50 * rsum * dt)

    elif scheme == CIRNumericalScheme.LOGNORMAL.value:

        x = np.exp(-a * dt)
        y = 1.0 - x

        for iPath in range(0, num_paths):

            r = r0
            rsum = r0

            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                mean = x * r + b * y
                var = sigma * sigma * y * (x * r + 0.50 * b * y) / a
                sig = np.sqrt(np.log(1.0 + var / (mean * mean)))
                r_prev = r
                r = mean * np.exp(-0.5 * sig * sig + sig * z[iStep - 1])
                rsum += (r + r_prev)

            zcb += np.exp(-0.5 * rsum * dt)

    elif scheme == CIRNumericalScheme.MILSTEIN.value:

        sigmasqrt_dt = sigma * np.sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0

        for iPath in range(0, num_paths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                r_prev = r
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrt_dt * np.sqrt(max(0.0, r))
                r = r + sigma2dt * (z[iStep - 1] ** 2 - 1.0)
                rsum += (r + r_prev)

            zcb += np.exp(-0.50 * rsum * dt)

    elif scheme == CIRNumericalScheme.KAHLJACKEL.value:

        bhat = b - sigma * sigma / 4.0 / a
        sqrt_dt = np.sqrt(dt)

        for iPath in range(0, num_paths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(num_steps - 1))

            for iStep in range(1, num_steps):
                beta = z[iStep - 1] / sqrt_dt
                rootr = np.sqrt(max(r, 1e-8))
                c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt / 4.0 / rootr
                r_prev = r
                r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt
                rsum += (r + r_prev)

            zcb += np.exp(-0.50 * rsum * dt)

    elif scheme == CIRNumericalScheme.EXACT.value:

        for iPath in range(0, num_paths):

            r = r0
            rsum = r

            for iStep in range(1, num_steps):
                r_prev = r
                r = draw(r, a, b, sigma, dt)
                rsum += (r + r_prev)

            zcb += np.exp(-0.50 * rsum * dt)

    zcb /= num_paths
    return zcb

###############################################################################
