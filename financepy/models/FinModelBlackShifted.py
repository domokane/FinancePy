##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from scipy.stats import norm

from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinGlobalTypes import FinOptionTypes

from ..finutils.FinMath import N

###############################################################################
# NOTE: Keeping this separate from SABR for the moment.
###############################################################################


class FinModelBlackShifted():
    ''' Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. This model also allows
    the distribution to be shifted to the negative in order to allow for
    negative interest rates. '''

    def __init__(self, volatility, shift, implementation=0):
        ''' Create FinModel black using parameters. '''
        self._volatility = volatility
        self._shift = shift
        self._implementation = 0
        self._numSteps = 0
        self._seed = 0
        self._param1 = 0
        self._param2 = 0

###############################################################################

    def value(self,
              forwardRate,   # Forward rate
              strikeRate,    # Strike Rate
              timeToExpiry,  # time to expiry in years
              df,            # Discount Factor to expiry date
              callOrPut):    # Call or put
        ''' Price a derivative using Black's model which values in the forward
        measure following a change of measure. The sign of the shift is the
        same as Matlab. '''

        s = self._shift
        f = forwardRate
        t = timeToExpiry
        k = strikeRate
        sqrtT = np.sqrt(t)
        vol = self._volatility

        d1 = np.log((f+s)/(k+s)) + vol * vol * t / 2
        d1 = d1 / (vol * sqrtT)
        d2 = d1 - vol * sqrtT

        if callOrPut == FinOptionTypes.EUROPEAN_CALL:
            return df * ((f+s) * N(d1) - (k+s) * N(d2))
        elif callOrPut == FinOptionTypes.EUROPEAN_PUT:
            return df * ((k+s) * N(-d2) - (f+s) * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VOLATILITY", self._volatility)
        s += labelToString("SHIFT", self._shift)
        s += labelToString("IMPLEMENTATION", self._implementation)
        s += labelToString("NUMSTEPS", self._numSteps)
        return s

###############################################################################
