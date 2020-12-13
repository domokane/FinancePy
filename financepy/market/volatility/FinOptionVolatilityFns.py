##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64

from ...finutils.FinMath import N

###############################################################################
# Parametric functions for option volatility to use in a Black-Scholes model
###############################################################################

@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def volFunctionClarke(params, F, K, texp):
    ''' Volatility Function in book by Iain Clarke generalised to allow for 
    higher than quadratic power. The '''

    x = np.log(F / K)
    sigma0 = np.exp(params[0])
    arg = x / (sigma0 * np.sqrt(texp))
    deltax = N(arg) - 0.50  # The -0.50 seems to be missing in book
    f = 0.0
    for i in range(0, len(params)):
        f += params[i] * (deltax ** i)

    return np.exp(f)

###############################################################################
