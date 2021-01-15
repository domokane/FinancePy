###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString

###############################################################################
# TODO: This should be deleted and replaced with FinEquityVolSurface


class FinEquityVolCurve():
    ''' Class to manage a smile or skew in volatility at a single maturity
    horizon. It fits the volatility using a polynomial. Includes analytics to
    extract the implied pdf of the underyling at maturity. THIS NEEDS TO BE
    SUBSTITUTED WITH FINEQUITYVOLSURFACE. '''

###############################################################################

    def __init__(self,
                 curveDate,
                 expiryDate,
                 strikes,
                 volatilities,
                 polynomial=3):

        if expiryDate <= curveDate:
            raise FinError("Expiry date before curve date.")

        if len(strikes) < 1:
            raise FinError("Volatility grid has zero length.")

        if testMonotonicity(strikes) is False:
            raise FinError("Strikes must be strictly monotonic.")

        numStrikes = len(strikes)
        numVols = len(volatilities)

        if numStrikes != numVols:
            raise FinError("Strike and volatility vectors not same length.")

        for i in range(1, numStrikes):
            if strikes[i] <= strikes[i - 1]:
                raise FinError("Grid Strikes are not in increasing order")

        self._curveDate = curveDate
        self._strikes = np.array(strikes)
        self._volatilities = np.array(volatilities)

        self._z = np.polyfit(self._strikes, self._volatilities, polynomial)
        self._f = np.poly1d(self._z)

###############################################################################

    def volatility(self, strike):
        ''' Return the volatility for a strike using a given polynomial
        interpolation. '''

        vol = self._f(strike)

        if vol.any() < 0.0:
            raise FinError("Negative volatility. Not permitted.")

        return vol

###############################################################################

    def calculatePDF():
        ''' calculate the probability density function of the underlying using
        the volatility smile or skew curve following the approach set out in
        Breedon and Litzenberger. '''
        pass

###############################################################################
