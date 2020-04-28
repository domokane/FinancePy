# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 15:28:41 2019

@author: Dominic O'Kane
"""
# TODO Fix this

import numpy as np
from scipy.stats import norm

from ...finutils.FinMath import N

###############################################################################
# NOTE: Need to convert option types to use enums.
# NOTE: Perhaps just turn this into a function rather than a class.
###############################################################################


class FinModelBlack():
    ''' Black's Model which prices call and put options in the forward
    measure according to the Black-Scholes equation. '''

    def value(self,
              forwardRate,
              strikeRate,
              timeToExpiry,
              sigma,
              callOrPut):
        ''' Price a derivative using Black's model which values in the forward
        measure following a change of measure. '''

        if strikeRate < 0:
            raise Exception("Error: Negative strike")

        if strikeRate == 0.0:
            strike = 1e-16

        sqrtT = np.sqrt(timeToExpiry)

        d1 = np.log(forwardRate / strike) + sigma * sigma * timeToExpiry / 2
        d1 = d1 / (sigma * sqrtT)
        d2 = d1 - sigma * sqrtT

        returnValue = 0.0

        if callOrPut == "C":
            returnValue = forwardRate * norm.cdf(d1) - strike * N(d2)
        elif callOrPut == "P":
            returnValue = strike * norm.cdf - forwardRate * N(-d1)
        else:
            raise Exception("Option type must be a Call(C) or a Put(P)")

        return returnValue

###############################################################################
###############################################################################
