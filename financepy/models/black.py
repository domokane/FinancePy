##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Greeks added thanks to Guillaume Lefieux
##############################################################################

# TODO Fix this

import numpy as np
from numba import njit, float64, float64

from ..utils.math import NVect, NPrimeVect
from ..utils.global_vars import gSmall
from ..utils.helpers import label_to_string
from ..utils.global_types import FinOptionTypes
from ..utils.error import FinError

###############################################################################
# TODO: Use Numba ?
###############################################################################

@njit(float64[:](float64, float64, float64, float64), fastmath=True, cache=True)
def calculateD1D2(f, t, k, v):

    t = np.maximum(t, gSmall)
    vol = np.maximum(v, gSmall)
    k = np.maximum(k, gSmall)
    sqrtT = np.sqrt(t)

    if f <= 0.0:
        raise FinError("Forward is zero.")

    if k <= 0.0:
        raise FinError("Strike is zero.")

    d1 = (np.log(f/k) + vol * vol * t / 2.0) / (vol * sqrtT)
    d2 = d1 - vol * sqrtT

    return np.array([d1, d2])

###############################################################################

class Black():
    """ Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. """

    def __init__(self, volatility, implementation=0):
        """ Create FinModel black using parameters. """
        self._volatility = volatility
        self._implementation = 0
        self._num_steps = 0
        self._seed = 0
        self._param1 = 0
        self._param2 = 0

###############################################################################


    def value(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # df RFR to expiry date
              call_or_put):    # Call or put
        """ Price a derivative using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        
        [d1, d2] = calculateD1D2(f, t, k, v)
        
        if call_or_put == FinOptionTypes.EUROPEAN_CALL:
            value = df * (f * NVect(d1) - k * NVect(d2))
        elif call_or_put == FinOptionTypes.EUROPEAN_PUT:
            value = df * (k * NVect(-d2) - f * NVect(-d1))
        else:
            raise FinError("Option type must be a European Call or Put")

        return value

###############################################################################


    def delta(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # RFR to expiry date
              call_or_put):    # Call or put
        """ Calculate delta using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility

        [d1, d2] = calculateD1D2(f, t, k, v)

        if call_or_put == FinOptionTypes.EUROPEAN_CALL:
            delta = df * NVect(d1)
        elif call_or_put == FinOptionTypes.EUROPEAN_PUT:
            delta = - df * NVect(-d1)
        else:
            raise FinError("Option type must be a European Call or Put")

        return delta

###############################################################################


    def gamma(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # RFR to expiry date
              call_or_put):    # Call or put
        """ Calculate gamma using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility

        [d1, d2] = calculateD1D2(f, t, k, v)

        sqrtT = np.sqrt(t)
        gamma = df * NPrimeVect(d1) / (f * v * sqrtT)
        return gamma

###############################################################################


    def theta(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # Discount Factor to expiry date
              call_or_put):    # Call or put
        """ Calculate theta using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        r = -np.log(df)/t

        [d1, d2] = calculateD1D2(f, t, k, v)

        sqrtT = np.sqrt(t)

        if call_or_put == FinOptionTypes.EUROPEAN_CALL:
            theta = df * (-(f * v * NPrimeVect(d1)) / (2*sqrtT) + r * f * NVect(d1)
                          - r * k * NVect(d2))        
        elif call_or_put == FinOptionTypes.EUROPEAN_PUT:
            theta = df * (-(f * v * NPrimeVect(d1)) / (2*sqrtT) - r * f * NVect(-d1)
                          + r * k * NVect(-d2))
        else:
            raise FinError("Option type must be a European Call or Put")

        return theta

###############################################################################

    def vega(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,  # df RFR to expiry date
              call_or_put):    # Call or put
        """ Price a derivative using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        v = self._volatility
        sqrtT = np.sqrt(t)
        
        [d1, d2] = calculateD1D2(f, t, k, v)
        
        if call_or_put == FinOptionTypes.EUROPEAN_CALL:
            vega = df * f * sqrtT * NPrimeVect(d1)
        elif call_or_put == FinOptionTypes.EUROPEAN_PUT:
            vega = df * f * sqrtT * NPrimeVect(d1)
        else:
            raise FinError("Option type must be a European Call or Put")

        return vega

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self._volatility)
        s += label_to_string("IMPLEMENTATION", self._implementation)
        s += label_to_string("NUMSTEPS", self._num_steps)
        return s

###############################################################################
