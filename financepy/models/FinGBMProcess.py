# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 15:05:58 2019

@author: Dominic
"""

import numpy as np
from math import sqrt, exp
from numba import njit, float64, int64
from ..finutils.FinHelperFunctions import labelToString

###############################################################################

@njit(float64[:, :](int64, int64, float64, float64, float64, float64, int64),
      cache=True, fastmath=True)
def getPaths(numPaths,
             numTimeSteps,
             t,
             mu,
             stockPrice,
             volatility,
             seed):

    np.random.seed(seed)
    dt = t / numTimeSteps
    vsqrtdt = volatility * sqrt(dt)
    m = exp((mu - volatility * volatility / 2.0) * dt)
    Sall = np.empty((2 * numPaths, numTimeSteps + 1))

    # This should be less memory intensive as we only generate randoms per step
    Sall[:, 0] = stockPrice
    for it in range(1, numTimeSteps + 1):
        g1D = np.random.standard_normal((numPaths))
        for ip in range(0, numPaths):
            w = np.exp(g1D[ip] * vsqrtdt)
            Sall[ip, it] = Sall[ip, it - 1] * m * w
            Sall[ip + numPaths, it] = Sall[ip + numPaths, it - 1] * m / w

    return Sall

##########################################################################


@njit(cache=True, fastmath=True)
def getPathsAssets(numAssets,
                   numPaths,
                   numTimeSteps,
                   t,
                   mus,
                   stockPrices,
                   volatilities,
                   betas,
                   seed):

    np.random.seed(seed)
    dt = t / numTimeSteps
    vsqrtdts = volatilities * sqrt(dt)
    m = np.exp((mus - volatilities * volatilities / 2.0) * dt)

    Sall = np.empty((2 * numPaths, numTimeSteps + 1, numAssets))

    g2D = np.random.standard_normal((numPaths, numTimeSteps + 1))
    g3D = np.random.standard_normal((numPaths, numTimeSteps + 1, numAssets))

    for ip in range(0, numPaths):
        for ia in range(0, numAssets):
            Sall[ip, 0, ia] = stockPrices[ia]
            Sall[ip + numPaths, 0, ia] = stockPrices[ia]

    for ip in range(0, numPaths):
        for it in range(1, numTimeSteps + 1):
            zmkt = g2D[ip, it]
            for ia in range(0, numAssets):
                zidio = g3D[ip, it, ia]
                z = betas[ia] * zmkt + sqrt(1.0 - betas[ia]**2) * zidio
                w = exp(z * vsqrtdts[ia])
                Sall[ip, it, ia] = Sall[ip, it - 1, ia] * m[ia] * w
                Sall[ip + numPaths, it, ia] = Sall[ip + \
                    numPaths, it - 1, ia] * m[ia] / w

    return Sall

##########################################################################


@njit(cache=True, fastmath=True)
def getAssets(numAssets,
              numPaths,
              t,
              mus,
              stockPrices,
              volatilities,
              betas,
              seed):

    np.random.seed(seed)
    vsqrtdts = volatilities * sqrt(t)
    m = np.exp((mus - volatilities * volatilities / 2.0) * t)

    Sall = np.empty((2 * numPaths, numAssets))

    g2D = np.random.standard_normal((numPaths))
    g3D = np.random.standard_normal((numPaths, numAssets))

    for ip in range(0, numPaths):
        zmkt = g2D[ip]
        for ia in range(0, numAssets):
            zidio = g3D[ip, ia]
            z = betas[ia] * zmkt + sqrt(1.0 - betas[ia]**2) * zidio
            w = exp(z * vsqrtdts[ia])
            Sall[ip, ia] = stockPrices[ia] * m[ia] * w
            Sall[ip + numPaths, ia] = stockPrices[ia] * m[ia] / w

    return Sall

##########################################################################


class FinGBMProcess():

    def getPaths(
            self,
            numPaths,
            numTimeSteps,
            t,
            mu,
            stockPrice,
            volatility,
            seed):
        paths = getPaths(numPaths, numTimeSteps,
                         t, mu, stockPrice, volatility, seed)
        return paths

    def getPathsAssets(self, numAssets, numPaths, numTimeSteps,
                       t, mus, stockPrices, volatilities, betas, seed):

        if numTimeSteps == 2:
            paths = getAssets(numAssets, numPaths,
                              t, mus, stockPrices, volatilities, betas, seed)
        else:
            paths = getPathsAssets(numAssets, numPaths, numTimeSteps,
                                   t, mus, stockPrices,
                                   volatilities, betas, seed)
        return paths

##########################################################################
