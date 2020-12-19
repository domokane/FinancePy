##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64

from ...finutils.FinMath import N


from enum import Enum

class FinVolFunctionTypes(Enum):
    CLARKE = 0
    SABR = 1


###############################################################################
# Parametric functions for option volatility to use in a Black-Scholes model
###############################################################################

@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def volFunctionClarke(params, f, k, t):
    ''' Volatility Function in book by Iain Clarke generalised to allow for 
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

@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def volFunctionSABR(params, f, k, t):
    ''' This is SABR but when beta = 1 '''
    
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
