# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np

##########################################################################
# TODO
# Inherit from FinCurve and add df method
# Put in a convention for the rate
# Use Frequency object
##########################################################################

#from ...market.curves.FinCurve import FinCurve
from ...finutils.FinInterpolate import FinInterpMethods

##########################################################################
##########################################################################


class FinPiecewiseLinearCurve():
    ''' Curve is made up of a series of sections assumed to each have a constant
    forward rate. This class needs to be checked carefully. '''

    def __init__(self, curveDate, times, values):
        ''' Curve is defined by a vector of increasing times and zero rates. '''

        if len(times) != len(values):
            raise ValueError("Times and rates vectors must have same length")

        if len(times) == 0:
            raise ValueError("Times and rates vectors must have length > 0")

        numTimes = len(times)
        for i in range(1, numTimes):
            if times[i - 1] >= times[i]:
                raise ValueError("Times must be strictly increasing")

        self._times = times
        self._values = values

##########################################################################

    def zero(self, t, interpolationMethod=FinInterpMethods.FLAT_FORWARDS):

        l_index = 0
        r_index = 0
        numTimes = len(self._times)
        for i in range(1, numTimes):
            if self._times[i] > t:
                l_index = i - 1
                r_index = i
                break

        interpolatedZero = 0.0

        if interpolationMethod == FinInterpMethods.FLAT:

            interpolatedZero = self._values[l_index]

        elif interpolationMethod == FinInterpMethods.LINEAR:

            t1 = self._times[l_index]
            t2 = self._times[r_index]
            r1 = self._values[l_index]
            r2 = self._values[r_index]
            interpolatedZero = r1 + (r2 - r1) * (t - t1) / (t2 - t1)

        elif interpolationMethod == FinInterpMethods.LOG:

            t1 = self._times[l_index]
            t2 = self._times[r_index]
            r1 = self._values[l_index]
            r2 = self._values[r_index]
            interpolatedZero = r1 * ((r2 / r1) ** (t - t1) / (t2 - t1))

        else:
            raise ValueError("Unknown interpolation scheme")

        return interpolatedZero

##########################################################################

    def fwd(self, t):
        # NEED TODO THIS
        fwdRate = self.r
        return fwdRate

##########################################################################

    def df(self,
           t,
           freq=0,  # This corresponds to continuous compounding
           interpolationMethod=FinInterpMethods.FLAT_FORWARDS):

        df = 1.0
        numTimes = len(self._times)
        index = 0

        for i in range(1, numTimes):
            if self._times[i] > t:
                index = i
                break

        for i in range(1, index):
            t = self._times[i]
            r = self._values[i]
            tau = self._times[i] - self._times[i - 1]

            if freq == 0:
                df = df * np.exp(-r * tau)
            elif freq > 0:
                df = df * ((1.0 + r / freq)**(-tau * freq))
            elif freq == -1:
                df = df * 1.0 / (1.0 + r * tau)

        tau = t - self._times[index]
        r = self.zero(t, interpolationMethod)

        if freq == 0:
            df = df * np.exp(-r * tau)
        elif freq > 0:
            df = df * ((1.0 + r / freq)**(-tau * freq))
        elif freq == -1:
            df = df * 1.0 / (1.0 + r * tau)

        return df

##########################################################################
