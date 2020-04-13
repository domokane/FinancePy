# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:10:12 2019

@author: Dominic
"""

from math import log, sqrt, exp
from numba import njit, float64, int64
import numpy as np

##########################################################################
# CIR Process
# dr = a(b-r) + sigma sqrt(r) dW
# Note that r can hit zero if 2.0 * a * b < sigma*sigma:
##########################################################################

# TO DO - DECIDE WHETHER TO OO MODEL
# CAN DO Z SCALING INSIDE NUMPY ?
# ANTITHETICS

from enum import Enum


class FinCIRNumericalScheme(Enum):
    EULER = 1
    LOGNORMAL = 2
    MILSTEIN = 3
    KAHLJACKEL = 4
    EXACT = 5  # SAMPLES EXACT DISTRIBUTION

##########################################################################

# THIS CLASS IS NOT USED BUT MAY BE USED IF WE CREATE AN OO FRAMEWORK


class FinModelRatesCIR():

    def __init__(self, a, b, sigma):
        self._a = a
        self._b = b
        self._sigma = sigma

##########################################################################


@njit(fastmath=True, cache=True)
def meanr(r0, a, b, t):
    ''' Mean value of a CIR process after time t '''
    mr = r0 * exp(-a * t) + b * (1.0 - exp(-a * t))
    return mr

##########################################################################


@njit(fastmath=True, cache=True)
def variancer(r0, a, b, sigma, t):
    ''' Variance of a CIR process after time t '''
    vr = r0 * sigma * sigma * (exp(-a * t) - exp(-2.0 * a * t)) / a
    vr += b * sigma * sigma * ((1.0 - exp(-a * t))**2) / 2.0 / a
    return vr

##########################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64),
    fastmath=True,
    cache=True)
def zeroPrice(r0, a, b, sigma, t):
    ''' Price of a zero coupon bond in CIR model. '''
    h = sqrt(a * a + 2.0 * sigma * sigma)
    denom = 2.0 * h + (a + h) * (exp(h * t) - 1.0)
    A = (2.0 * h * exp((a + h) * t / 2.0) /
         denom)**(2.0 * a * b / sigma / sigma)
    B = 2.0 * (exp(h * t) - 1.0) / denom
    zcb = A * exp(-r0 * B)
    return zcb

##########################################################################


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
    ''' Draw a next rate from the CIR model in Monte Carlo. '''
    sigma2 = sigma * sigma
    d = 4.0 * a * b / sigma2
    l = 4.0 * a * exp(-a * dt) / sigma2 / (1.0 - exp(-a * dt)) * rt
    c = sigma2 * (1.0 - exp(-a * dt)) / 4.0 / a

    if d > 1:
        Z = np.random.normal()
        X = np.random.chisquare(d - 1)
        r = c * (X + (Z + sqrt(l))**2)
    else:
        N = np.random.poisson(l / 2.0)
        X = np.random.chisquare(d + 2 * N)
        r = c * X

    return r

##########################################################################


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
def ratePath_MC(r0, a, b, sigma, t, dt, seed, scheme):
    ''' Generate a path of CIR rates using a number of numerical schemes. '''

    np.random.seed(seed)
    numSteps = int(t / dt)
    ratePath = np.zeros(numSteps)
    ratePath[0] = r0
    numPaths = 1

    if scheme == FinCIRNumericalScheme.EULER.value:

        sigmasqrtdt = sigma * sqrt(dt)

        for iPath in range(0, numPaths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrtdt * sqrt(max(r, 0.0))
                ratePath[iStep] = r

    elif scheme == FinCIRNumericalScheme.LOGNORMAL.value:

        x = exp(-a * dt)
        y = 1.0 - x

        for iPath in range(0, numPaths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                mean = x * r + b * y
                var = sigma * sigma * y * (x * r + 0.50 * b * y) / a
                sig = sqrt(log(1.0 + var / (mean * mean)))
                r = mean * exp(-0.5 * sig * sig + sig * z[iStep - 1])
                ratePath[iStep] = r

    elif scheme == FinCIRNumericalScheme.MILSTEIN.value:

        sigmasqrtdt = sigma * sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0

        for iPath in range(0, numPaths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrtdt * sqrt(max(0.0, r))
                r = r + sigma2dt * (z[iStep - 1]**2 - 1.0)
                ratePath[iStep] = r

    elif scheme == FinCIRNumericalScheme.KAHLJACKEL.value:

        bhat = b - sigma * sigma / 4.0 / a
        sqrtdt = sqrt(dt)

        for iPath in range(0, numPaths):

            r = r0
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                beta = z[iStep - 1] / sqrtdt
                rootr = sqrt(max(r, 1e-8))
                c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt / 4.0 / rootr
                r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt
                ratePath[iStep] = r

    elif scheme == FinCIRNumericalScheme.EXACT.value:

        for iPath in range(0, numPaths):

            r = r0

            for iStep in range(1, numSteps):
                r = draw(r, a, b, sigma, dt)
                ratePath[iStep] = r

    return ratePath

##########################################################################


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
def zeroPrice_MC(r0, a, b, sigma, t, dt, numPaths, seed, scheme):
    ' Determine the CIR zero price using Monte Carlo. '''

    if t == 0.0:
        return 1.0

    np.random.seed(seed)

    numSteps = int(t / dt)
    zcb = 0.0

    if scheme == FinCIRNumericalScheme.EULER.value:

        sigmasqrtdt = sigma * sqrt(dt)

        for iPath in range(0, numPaths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                r_prev = r
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrtdt * sqrt(max(r, 0.0))
                rsum += (r + r_prev)

            zcb += exp(-0.50 * rsum * dt)

    elif scheme == FinCIRNumericalScheme.LOGNORMAL.value:

        x = exp(-a * dt)
        y = 1.0 - x

        for iPath in range(0, numPaths):

            r = r0
            rsum = r0

            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):

                mean = x * r + b * y
                var = sigma * sigma * y * (x * r + 0.50 * b * y) / a
                sig = sqrt(log(1.0 + var / (mean * mean)))
                r_prev = r
                r = mean * exp(-0.5 * sig * sig + sig * z[iStep - 1])
                rsum += (r + r_prev)

            zcb += exp(-0.5 * rsum * dt)

    elif scheme == FinCIRNumericalScheme.MILSTEIN.value:

        sigmasqrtdt = sigma * sqrt(dt)
        sigma2dt = sigma * sigma * dt / 4.0

        for iPath in range(0, numPaths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                r_prev = r
                r = r + a * (b - r) * dt + \
                    z[iStep - 1] * sigmasqrtdt * sqrt(max(0.0, r))
                r = r + sigma2dt * (z[iStep - 1]**2 - 1.0)
                rsum += (r + r_prev)

            zcb += exp(-0.50 * rsum * dt)

    elif scheme == FinCIRNumericalScheme.KAHLJACKEL.value:

        bhat = b - sigma * sigma / 4.0 / a
        sqrtdt = sqrt(dt)

        for iPath in range(0, numPaths):

            r = r0
            rsum = r
            z = np.random.normal(0.0, 1.0, size=(numSteps - 1))

            for iStep in range(1, numSteps):
                beta = z[iStep - 1] / sqrtdt
                rootr = sqrt(max(r, 1e-8))
                c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt / 4.0 / rootr
                r_prev = r
                r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt
                rsum += (r + r_prev)

            zcb += exp(-0.50 * rsum * dt)

    elif scheme == FinCIRNumericalScheme.EXACT.value:

        for iPath in range(0, numPaths):

            r = r0
            rsum = r

            for iStep in range(1, numSteps):
                r_prev = r
                r = draw(r, a, b, sigma, dt)
                rsum += (r + r_prev)

            zcb += exp(-0.50 * rsum * dt)

    zcb /= numPaths
    return zcb

##########################################################################
