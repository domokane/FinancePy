# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinInterpolate import interpolate

###############################################################################
###############################################################################


class FinCurve():
    ''' Class to manage curves from which other curve class inherit. '''

##########################################################################

    def __init__(self, times, values):
        ''' Create curve as a vector of times and values of same length. '''

        if len(times) < 1:
            print("Times has zero length")
            return

        if len(times) != len(values):
            print("Times and Values are not the same")
            return

        num_times = len(times)

        if times[0] < 0.0:
            print("First time is negative")
            return

        for i in range(1, num_times):

            if times[i] <= times[i - 1]:
                print("Times are not sorted in increasing order")

        self._times = times
        self._values = values

##########################################################################

    def value(self, t):
        ''' get the value (in this context it is a discount factor) by
        interpolating according to a specified interpolation method. '''

        v = interpolate(t, self._times, self._values, self._method)
        return v

##########################################################################
