##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from scipy.stats import norm

from ..finutils.FinHelperFunctions import labelToString
from ..finutils.FinGlobalTypes import FinOptionTypes

###############################################################################
# NOTE: Need to convert option types to use enums.
# NOTE: Perhaps just turn this into a function rather than a class.
###############################################################################


class FinModelBachelier():
    ''' Bachelier's Model which prices call and put options in the forward
    measure assuming the underlying rate follows a normal process. '''

    def __init__(self, volatility):
        ''' Create FinModel black using parameters. '''
        self._volatility = volatility

###############################################################################

    def value(self,
              forwardRate,   # Forward rate F
              strikeRate,    # Strike Rate K
              timeToExpiry,  # Time to Expiry (years)
              df,            # Discount Factor to expiry date
              callOrPut):    # Call or put
        ''' Price a call or put option using Bachelier's model. '''

        f = forwardRate
        t = timeToExpiry
        k = strikeRate
        sqrtT = np.sqrt(t)
        vol = self._volatility
        d = (f-k) / (vol * sqrtT)

        if callOrPut == FinOptionTypes.EUROPEAN_CALL:
            return df * ((f - k) * norm.cdf(d) + vol * sqrtT * norm.pdf(d))
        elif callOrPut == FinOptionTypes.EUROPEAN_PUT:
            return df * ((k - f) * norm.cdf(-d) + vol * sqrtT * norm.pdf(d))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

        return 999

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VOLATILITY", self._volatility)
        return s

###############################################################################
