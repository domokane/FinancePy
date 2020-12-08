##############################################################################
# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi
##############################################################################

import numpy as np
from scipy import optimize
from numba import njit, float64, int64, vectorize

from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinGlobalVariables import gSmall
from ..finutils.FinMath import NVect, NPrimeVect, N
from ..finutils.FinError import FinError

###############################################################################
# Analytical Black Scholes model implementation and approximations
###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsValue(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif optionTypeValue == FinOptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    vsqrtT = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    d2 = d1 - vsqrtT
    value = phi * ss * NVect(phi*d1)
    value += - phi * kk * NVect(phi*d2)
    return value

###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsDelta(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif optionTypeValue == FinOptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    vsqrtT = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    delta = phi * np.exp(-q*t) * NVect(phi*d1)
    return delta

###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsGamma(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    vsqrtT = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    gamma = np.exp(-q*t) * NPrimeVect(d1) / s / vsqrtT
    return gamma

###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsVega(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    sqrtT = np.sqrt(t)
    vsqrtT = v * sqrtT
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    vega = ss * sqrtT * NPrimeVect(d1)
    return vega

###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsTheta(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif optionTypeValue == FinOptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    sqrtT = np.sqrt(t)
    vsqrtT = v * sqrtT
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    d2 = d1 - vsqrtT
    theta = - ss * NPrimeVect(d1) * v / 2.0 / sqrtT
    theta = theta - phi * r * k * np.exp(-r*t) * NVect(phi*d2)
    theta = theta + phi * q * ss * NVect(phi*d1)
    return theta

###############################################################################

@vectorize([float64(float64, float64, float64, float64, float64, float64, int64)],
           fastmath=True, cache=True)
def bsRho(s, t, k, r, q, v, optionTypeValue):
    ''' Price a derivative using Black-Scholes model. ''' 

    if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif optionTypeValue == FinOptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    sqrtT = np.sqrt(t)
    vsqrtT = v * sqrtT
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / vsqrtT + vsqrtT / 2.0
    d2 = d1 - vsqrtT
    rho = phi * k * t * np.exp(-r*t) * NVect(phi*d2)
    return rho

###############################################################################
###############################################################################
# This module contains a number of analytical approximations for the price of
# an American style option starting with Barone-Adesi-Whaley
# https://deriscope.com/docs/Barone_Adesi_Whaley_1987.pdf
###############################################################################
###############################################################################

def _fcall(si, *args):
    ''' Function to determine ststar for pricing American call options. '''

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    M = 2.0 * r / v2
    W = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q2 = (1.0 - W + np.sqrt((W - 1.0)**2 + 4.0 * M/K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, +1) 
    objFn = objFn - (1.0 - np.exp(-q*t) * N(d1)) * si / q2
    return objFn

###############################################################################

def _fput(si, *args):
    ''' Function to determine sstar for pricing American put options. '''

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    W = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q1 = (1.0 - W - np.sqrt((W - 1.0)**2 + 4.0 * K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, -1)
    objFn = objFn - (1.0 - np.exp(-q*t) * N(-d1)) * si / q1
    return objFn

###############################################################################

def bawValue(stockPrice: float, 
             timeToExpiry: float, 
             strikePrice: float, 
             riskFreeRate: float, 
             dividendRate: float, 
             volatility: float, 
             phi:int):
    ''' American Option Pricing Approximation using the Barone-Adesi-Whaley
    approximation for the Black Scholes Model '''

    s = stockPrice
    k = strikePrice
    t = timeToExpiry
    r = riskFreeRate
    q = dividendRate
    v = volatility
    b = riskFreeRate - dividendRate

    if phi == 1:

        if b >= r:
            return bsValue(s, t, k, r, q, v, +1)

        argtuple = (t, k, r, q, v)

        sstar = optimize.newton(_fcall, x0=s, fprime=None, args=argtuple,
                                tol=1e-7, maxiter=50, fprime2=None)
                    
        M = 2.0 * r / (v*v)
        W = 2.0 * b / (v*v) 
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar/k) + (b + v*v/ 2.0) * t) / (v * np.sqrt(t))
        q2 = (-1.0 * (W - 1.0) + np.sqrt((W - 1.0)**2 + 4.0 * M/K)) / 2.0
        A2 = (sstar / q2) * (1.0 - np.exp(-q * t) * N(d1))

        if s < sstar:
            return bsValue(s, t, k, r, q, v, +1) + A2 * ((s/sstar)**q2)
        else:
            return s - k

    elif phi == -1:

        argtuple = (t, k, r, q, v)
        sstar = optimize.newton(_fput, x0=s, fprime=None, args=argtuple,
                                tol=1e-7, maxiter=50, fprime2=None)
        v2 = v * v
        
        M = 2.0 * r / v2
        W = 2.0 * b / v2
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (-1.0 * (W - 1.0) - np.sqrt((W - 1.0)**2 + 4.0 * M/K)) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * N(-d1))

        if s > sstar:
            return bsValue(s, t, k, r, q, v, -1) + a1 * ((s/sstar)**q1)
        else:
            return k - s

    else:
        
        raise FinError("Phi must equal 1 or -1.")

###############################################################################

if __name__ == '__main__':
    # spot_price, strikePrice, timeToExpiry, r, b, vol, phi

    # Checking against table 3-1 in Haug 
    k = 100.0
    r = 0.10
    q = 0.10
    
    for t in [0.1, 0.5]:
        for v in [0.15, 0.25, 0.35]:
            for s in [90.0, 100.0, 110.0]:
                bawPrice = bawValue(s, t, k, r, q, v, +1)
                print("%9.5f %9.5f %9.5f %9.5f"% (s, t, v, bawPrice))

