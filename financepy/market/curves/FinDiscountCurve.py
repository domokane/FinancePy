# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from math import log
import numpy as np

from ...finutils.FinHelperFunctions import inputTime, inputFrequency
from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinMath import testMonotonicity
from .FinInterpolate import interpolate

##########################################################################


class FinDiscountCurve():

##########################################################################

    def __init__(self, curveDate, times, values, interpMethod):

        # Validate curve
        if len(times) < 1:
            raise FinError("Times has zero length")

        if len(times) != len(values):
            raise FinError("Times and Values are not the same")

        if times[0] != 0.0:
            raise FinError("First time is not zero")

        if values[0] != 1.0:
            raise FinError("First value is not zero" + str(values[0]))

        if testMonotonicity(times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._curveDate = curveDate
        self._times = np.array(times)
        self._values = np.array(values)
        self._interpMethod = interpMethod

##########################################################################

    def zeroRate(self, dt, compoundingFreq=-1):
        ''' Calculate the zero rate to maturity date. '''
        t = inputTime(dt, self)
        f = inputFrequency(compoundingFreq)
        df = self.df(t)

        if f == 0:  # Simple interest
            zeroRate = (1.0/df-1.0)/t
        if f == -1:  # Continuous
            zeroRate = -np.log(df) / t
        else:
            zeroRate = (df**(-1.0/t/f) - 1) * f
        return zeroRate

##########################################################################

    def df(self, dt):
        t = inputTime(dt, self)
        z = interpolate(t, self._times, self._values, self._interpMethod.value)
        return z

##########################################################################

    def survProb(self, dt):
        t = inputTime(dt, self)
        q = interpolate(t, self._times, self._values, self._interpMethod.value)
        return q

##########################################################################

    def fwd(self, dt):
        ''' Calculate the continuous forward rate at the forward date. '''
        t = inputTime(dt, self)
        dt = 0.000001
        df1 = self.df(t)
        df2 = self.df(t+dt)
        fwd = np.log(df1/df2)/dt
        return fwd

##########################################################################

    def fwdRate(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according to the specified
        day count convention. '''

        if date1 < self._curveDate:
            raise ValueError("Date1 before curve value date.")

        if date2 < date1:
            raise ValueError("Date2 must not be before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1, date2)
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / yearFrac
        return fwd

##########################################################################

    def print(self):
        numPoints = len(self._times)
        print("TIMES,DISCOUNT FACTORS")
        for i in range(0, numPoints):
            print("%10.7f,%10.7f" % (self._times[i], self._values[i]))

#######################################################################
