##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64

from ..finutils.FinMath import N
from ..finutils.FinError import FinError

###############################################################################
# Parametric functions for option volatility to use in a Black-Scholes model
###############################################################################

from enum import Enum

class FinVolFunctionTypes(Enum):
    CLARK = 0
    SABR3 = 1
    BBG = 2
    SABR = 3
    CLARK5= 4

###############################################################################


@njit(float64(float64[:], float64, float64, float64),
      fastmath=True, cache=True)
def volFunctionClark(params, f, k, t):
    ''' Volatility Function in book by Iain Clark generalised to allow for 
    higher than quadratic power. Care needs to be taken to avoid overfitting. 
    The exact reference is Clarke Page 59. '''

    if f < 0.0:
        print("f:", f)
        raise FinError("Forward is negative")

    if k < 0.0:
        print("k:", k)
        raise FinError("Strike is negative")

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


