##############################################################################
# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi
##############################################################################

import numpy as np
from financepy.finutils.FinMath import normcdf_fast, normpdf

from financepy.models.FinModelBlackScholes import bsValue
from scipy import optimize

normcdf = normcdf_fast
NP = normpdf

###############################################################################
# This module contains a number of analytical approximations for the price of
# an American style option starting with Barone-Adesi-Whaley
# https://deriscope.com/docs/Barone_Adesi_Whaley_1987.pdf
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
    N = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q2 = (1.0 - N + np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, +1) 
    objFn = objFn - (1.0 - np.exp(-q*t) * normcdf(d1)) * si / q2
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

    N = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q1 = (1.0 - N - np.sqrt((N - 1.0)**2 + 4.0 * K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, -1)
    objFn = objFn - (1.0 - np.exp(-q*t) * normcdf(-d1)) * si / q1
    return objFn

###############################################################################

def americanOptionBAW(self, 
                      stockPrice: float, 
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
    b = self._riskFreeRate - self._dividendRate

    if self._phi == 1:

        if b >= r:
            return bsValue(s, t, k, r, q, v, +1)

        argtuple = (t, k, r, q, v)

        sstar = optimize.newton(_fcall, x0=s, fprime=None, args=argtuple,
                                tol=1e-7, maxiter=50, fprime2=None)
                    
        M = 2.0 * r / (v*v)
        N = 2.0 * b / (v*v) 
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar/k) + (b + v*v/ 2.0) * t) / (v * np.sqrt(t))
        q2 = (-1.0 * (N - 1.0) + np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
        A2 = (sstar / q2) * (1.0 - np.exp(-q * t) * normcdf(d1))

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
        N = 2.0 * b / v2
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (-1.0 * (N - 1.0) - np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * normcdf(-d1))

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
                bawPrice = valueAmericanOptionBAW(s, t, k, r, q, v, +1)
                print("%9.5f %9.5f %9.5f %9.5f"% (s, t, v, bawPrice))

