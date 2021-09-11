##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from scipy.stats import norm

from ..utils.helpers import label_to_string
from ..utils.global_types import OptionTypes

###############################################################################
# NOTE: Need to convert option types to use enums.
# NOTE: Perhaps just turn this into a function rather than a class.
###############################################################################


class Bachelier():
    """ Bachelier's Model which prices call and put options in the forward
    measure assuming the underlying rate follows a normal process. """

    def __init__(self, volatility):
        """ Create FinModel black using parameters. """
        self._volatility = volatility

###############################################################################

    def value(self,
              forward_rate,   # Forward rate F
              strike_rate,    # Strike Rate K
              time_to_expiry,  # Time to Expiry (years)
              df,            # Discount Factor to expiry date
              call_or_put):    # Call or put
        """ Price a call or put option using Bachelier's model. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        sqrtT = np.sqrt(t)
        vol = self._volatility
        d = (f-k) / (vol * sqrtT)

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * ((f - k) * norm.cdf(d) + vol * sqrtT * norm.pdf(d))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * ((k - f) * norm.cdf(-d) + vol * sqrtT * norm.pdf(d))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

        return 999

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self._volatility)
        return s

###############################################################################
