##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64

from ...finutils.FinMath import N

###############################################################################

from enum import Enum

class FinVolFunctionTypes(Enum):
    CLARK = 0
    SABR = 1
    BBG = 2

###############################################################################
# Parametric functions for option volatility to use in a Black-Scholes model
###############################################################################


@njit(float64(float64[:], float64, float64, float64),
      fastmath=True, cache=True)
def volFunctionClark(params, f, k, t):
    ''' Volatility Function in book by Iain Clark generalised to allow for 
    higher than quadratic power. Care needs to be taken to avoid overfitting. 
    The exact reference is Clarke Page 59. '''

    x = np.log(f/k)
    sigma0 = np.exp(params[0])
    arg = x / (sigma0 * np.sqrt(t))
    deltax = N(arg) - 0.50  # The -0.50 seems to be missing in book
    f = 0.0
    for i in range(0, len(params)):
        f += params[i] * (deltax ** i)

    return np.exp(f)

###############################################################################


@njit(float64(float64[:], float64, float64, float64), 
      fastmath=True, cache=True)
def volFunctionBloomberg(params, f, k, t):
    ''' Volatility Function similar to the one used by Bloomberg. It is 
    a quadratic function in the spot delta of the option. It can therefore 
    go negative so it requires a good initial guess when performing the 
    fitting to avoid this happening. The first parameter is the quadratic 
    coefficient i.e. sigma(K) = a * D * D + b * D + c where a = params[0], 
    b = params[1], c = params[2] and D is the spot delta.'''
 
    numParams = len(params)

    # Rather than pass in the ATM vol, I imply it from the delta=0.50 curve
    sigma = 0.0
    for i in range(0, len(params)):
        pwr = numParams - i - 1
        sigma += params[i] * ((0.50) ** pwr)

    vsqrtt = sigma * np.sqrt(t)
    
    d1 = np.log(f/k)/ vsqrtt + vsqrtt/2.0
    delta = N(d1)

    v = 0.0
    for i in range(0, len(params)):
        pwr = numParams - i - 1
        v += params[i] * (delta ** pwr)

    return v

###############################################################################


@njit(float64(float64[:], float64, float64, float64), 
      fastmath=True, cache=True)
def volFunctionSABR(params, f, k, t):
    ''' This is the SABR function with the exponent beta set equal to 1. The 
    first parameter is alpha, then nu and the third parameter is rho. '''
    
    alpha = params[0]
    nu = params[1]
    rho = params[2]

    if rho > 1.0:
        rho = 0.99
    
    if rho < -1.0:
        rho = -0.99

    m = f / k

    if abs(m - 1.0) > 1e-6:

        sigma = 1.0
        numTerm1 = 0.0
        numTerm2 = rho * nu * alpha / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        logM = np.log(m)
        z = nu / alpha * logM
        denom = 1.0
        x = np.log((np.sqrt(1.0 - 2.0*rho*z + z**2.0) + z - rho)/(1.0 - rho))
        sigma = num*z/(denom*x)

    else:
        # when the option is at the money
        numTerm1 = 0.0
        numTerm2 = rho * nu * alpha / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        denom = 1.0
        sigma = num / denom

    return sigma

###############################################################################
