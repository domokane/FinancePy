##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from scipy.stats import norm

from ..finutils.FinMath import N
from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinOptionTypes import FinOptionTypes

###############################################################################
# NOTE: Need to convert option types to use enums.
# NOTE: Perhaps just turn this into a function rather than a class.
###############################################################################


class FinModelBlack():
    ''' Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. '''

    def __init__(self, volatility, implementation=0):
        ''' Create FinModel black using parameters. '''
        self._volatility = volatility
        self._implementation = 0
        self._numSteps = 0
        self._seed = 0
        self._param1 = 0
        self._param2 = 0

###############################################################################

    def value(self,
              forwardRate,   # Forward rate F
              strikeRate,    # Strike Rate K
              timeToExpiry,  # Time to Expiry (years)
              df,            # Discount Factor to expiry date
              callOrPut):    # Call or put
        ''' Price a derivative using Black's model which values in the forward
        measure following a change of measure. '''

        f = forwardRate
        t = timeToExpiry
        k = strikeRate
        sqrtT = np.sqrt(t)
        vol = self._volatility

        d1 = (np.log((f)/(k)) + vol * vol * t / 2) / (vol * sqrtT)
        d2 = d1 - vol*sqrtT

        if callOrPut == FinOptionTypes.EUROPEAN_CALL:
            return df * (f * norm.cdf(d1) - k * N(d2))
        elif callOrPut == FinOptionTypes.EUROPEAN_PUT:
            return df * (k * norm.cdf(-d2) - f * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

        return 999

###############################################################################

    def __repr__(self):
        s = "FINMODELBLACK"
        s += labelToString("VOLATILITY", self._volatility)
        s += labelToString("IMPLEMENTATION", self._implementation)
        s += labelToString("NUMSTEPS", self._numSteps)
        return s

###############################################################################
