##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Greeks added thanks to Guillaume Lefieux
##############################################################################

# TODO Fix this

import numpy as np
from numba import njit, float64

from ..utils.math import n_vect, n_prime_vect
from ..utils.global_vars import gSmall
from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes
from ..utils.error import FinError

###############################################################################
# TODO: Use Numba ?
###############################################################################


@njit(float64[:](float64, float64, float64, float64), fastmath=True, cache=True)
def calculate_d1_d2(f, t, k, v):

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

        [d1, d2] = calculate_d1_d2(f, t, k, v)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            value = df * (f * n_vect(d1) - k * n_vect(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            value = df * (k * n_vect(-d2) - f * n_vect(-d1))
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

        [d1, d2] = calculate_d1_d2(f, t, k, v)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            delta = df * n_vect(d1)
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            delta = - df * n_vect(-d1)
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

        [d1, d2] = calculate_d1_d2(f, t, k, v)

        sqrtT = np.sqrt(t)
        gamma = df * n_prime_vect(d1) / (f * v * sqrtT)
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

        [d1, d2] = calculate_d1_d2(f, t, k, v)

        sqrtT = np.sqrt(t)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            theta = df * (-(f * v * n_prime_vect(d1)) / (2 * sqrtT) + r * f * n_vect(d1)
                          - r * k * n_vect(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            theta = df * (-(f * v * n_prime_vect(d1)) / (2 * sqrtT) - r * f * n_vect(-d1)
                          + r * k * n_vect(-d2))
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

        [d1, d2] = calculate_d1_d2(f, t, k, v)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            vega = df * f * sqrtT * n_prime_vect(d1)
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            vega = df * f * sqrtT * n_prime_vect(d1)
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
