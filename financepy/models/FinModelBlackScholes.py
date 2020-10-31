##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from numba import njit

from scipy.stats import norm
from ..finutils.FinGlobalVariables import gSmall
N = norm.cdf

###############################################################################
# This is intended to be a fast calculator and validation is left to calling
# functions.
###############################################################################


def bsValue(s:float, 
            t:float, 
            k:float, 
            r:float, 
            q:float, 
            v:float, 
            phi:int):   # +1 for call, -1 for put
    ''' Price a derivative using Black-Scholes model where 
    phi is +1 for a call, and
    phi is -1 for a put.'''

    k = np.maximum(k, gSmall)
    t = np.maximum(t, gSmall)
    v = np.maximum(v, gSmall)

    sqrtT = np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v / sqrtT + v * sqrtT / 2.0
    d2 = d1 - v * sqrtT
    v = phi * ss * N(phi*d1) - phi * kk * N(phi*d2)
    return v

###############################################################################
