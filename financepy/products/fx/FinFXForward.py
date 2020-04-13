# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from scipy import optimize
from scipy.stats import norm

#from ...finutils.FinMath import N, nprime
from ...finutils.FinDate import FinDate
from ...finutils.FinMath import nprime
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError

N = norm.cdf

###############################################################################


class FinFXForward():

    def __init__(self,
                 expiryDate,
                 forwardFXRate): # value of a unit of foreign in domestic currency

        self._expiryDate = expiryDate
        self._fowardFXRate = forwardFXRate

###############################################################################

    def value(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        forDf = forDiscountCurve.df(t)

        v = spotFXRate * forDF / domDF
        return v

###############################################################################

    def forward(self,
              valueDate,
              spotFXRate, # value of a unit of foreign in domestic currency
              domDiscountCurve,
              forDiscountCurve):

        if type(valueDate) == FinDate:
            t = (self._expiryDate - valueDate) / gDaysInYear
        else:
            t = valueDate

        if np.any(spotFXRate <= 0.0):
            raise FinError("spotFXRate must be greater than zero.")

        if np.any(t < 0.0):
            raise FinError("Time to expiry must be positive.")

        t = np.maximum(t, 1e-10)

        domDf = domDiscountCurve.df(t)
        forDf = forDiscountCurve.df(t)

        fwd = spotFXRate * forDF / domDF
        return fwd

###############################################################################
