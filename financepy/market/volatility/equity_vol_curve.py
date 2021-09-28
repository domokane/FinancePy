###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.math import test_monotonicity

###############################################################################
# TODO: This should be deleted and replaced with equity_vol_surface


class EquityVolCurve():
    """ Class to manage a smile or skew in volatility at a single maturity
    horizon. It fits the volatility using a polynomial. Includes analytics to
    extract the implied pdf of the underyling at maturity. THIS NEEDS TO BE
    SUBSTITUTED WITH FINEQUITYVOLSURFACE. """

###############################################################################

    def __init__(self,
                 curve_date,
                 expiry_date,
                 strikes,
                 volatilities,
                 polynomial=3):

        if expiry_date <= curve_date:
            raise FinError("Expiry date before curve date.")

        if len(strikes) < 1:
            raise FinError("Volatility grid has zero length.")

        if test_monotonicity(strikes) is False:
            raise FinError("Strikes must be strictly monotonic.")

        num_strikes = len(strikes)
        num_vols = len(volatilities)

        if num_strikes != num_vols:
            raise FinError("Strike and volatility vectors not same length.")

        for i in range(1, num_strikes):
            if strikes[i] <= strikes[i - 1]:
                raise FinError("Grid Strikes are not in increasing order")

        self._curve_date = curve_date
        self._strikes = np.array(strikes)
        self._volatilities = np.array(volatilities)

        self._z = np.polyfit(self._strikes, self._volatilities, polynomial)
        self._f = np.poly1d(self._z)

###############################################################################

    def volatility(self, strike):
        """ Return the volatility for a strike using a given polynomial
        interpolation. """

        vol = self._f(strike)

        if vol.any() < 0.0:
            raise FinError("Negative volatility. Not permitted.")

        return vol

###############################################################################

    def calculate_pdf():
        """ calculate the probability density function of the underlying using
        the volatility smile or skew curve following the approach set out in
        Breedon and Litzenberger. """
        pass

###############################################################################
