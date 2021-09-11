##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np

from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes

from ..utils.math import N

###############################################################################
# NOTE: Keeping this separate from SABR for the moment.
###############################################################################


class BlackShifted():
    """ Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. This model also allows
    the distribution to be shifted to the negative in order to allow for
    negative interest rates. """

    def __init__(self, volatility, shift, implementation=0):
        """ Create FinModel black using parameters. """
        self._volatility = volatility
        self._shift = shift
        self._implementation = 0
        self._num_steps = 0
        self._seed = 0
        self._param1 = 0
        self._param2 = 0

###############################################################################

    def value(self,
              forward_rate,   # Forward rate
              strike_rate,    # Strike Rate
              time_to_expiry,  # time to expiry in years
              df,            # Discount Factor to expiry date
              call_or_put):    # Call or put
        """ Price a derivative using Black's model which values in the forward
        measure following a change of measure. The sign of the shift is the
        same as Matlab. """

        s = self._shift
        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        sqrtT = np.sqrt(t)
        vol = self._volatility

        d1 = np.log((f+s)/(k+s)) + vol * vol * t / 2
        d1 = d1 / (vol * sqrtT)
        d2 = d1 - vol * sqrtT

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * ((f+s) * N(d1) - (k + s) * N(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * ((k+s) * N(-d2) - (f + s) * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self._volatility)
        s += label_to_string("SHIFT", self._shift)
        s += label_to_string("IMPLEMENTATION", self._implementation)
        s += label_to_string("NUMSTEPS", self._num_steps)
        return s

###############################################################################
