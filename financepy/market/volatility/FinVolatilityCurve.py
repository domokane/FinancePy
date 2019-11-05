# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np

from ...finutils.FinInterpolate import FinInterpMethods, interpolate
from ...finutils.FinError import FinError

##########################################################################
# TODO: This is unfinished


class FinVolatilityCurve():

    def __init__(self,
                 curveDate,
                 strikeVector,
                 volatilityVector):

        if len(strikeVector) < 1:
            raise FinError("Volatility grid has zero length")

        numStrikes = len(strikeVector)
        numVols = len(volatilityVector)

        if numStrikes != numVols:
            raise FinError("Strike and volatility vectors not same length.")

        for i in range(1, numStrikes):
            if strikeVector[i] <= strikeVector[i - 1]:
                raise FinError("Grid Strikes are not in increasing order")

        self._curveDate = curveDate
        self._strikeVector = np.array(strikeVector)
        self._volatilityVector = np.array(volatilityVector)

###############################################################################

    def volatility(self,
                   strike,
                   interpMethod=FinInterpMethods.PIECEWISE_LINEAR):
        ''' Return the volatility for a strike using a given interpolation. '''

        vol = interpolate(strike,
                          self._strikeVector,
                          self._volatilityVector,
                          interpMethod.value)

        if vol < 0.0:
            raise ValueError("Negative volatility")

        return vol

###############################################################################

    def fitPolynomial(self,
                      strike,
                      interpMethod=FinInterpMethods.PIECEWISE_LINEAR):

        vol = interpolate(strike,
                          self._strikeVector,
                          self._volsVector,
                          interpMethod.value)

        return vol
