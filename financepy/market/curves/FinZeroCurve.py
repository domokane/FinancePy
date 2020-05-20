# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np

from ...finutils.FinHelperFunctions import inputTime, inputFrequency

from ...finutils.FinFrequency import FinFrequency, FinFrequencyTypes
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes
from ...finutils.FinMath import testMonotonicity
from .FinInterpolate import interpolate, FinInterpMethods
from ...finutils.FinHelperFunctions import labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve


###############################################################################
# TODO: Fix up __repr__ function
###############################################################################


class FinZeroCurve():
    ''' This is a curve calculated from a set of times and zero rates.
    '''

###############################################################################

    def __init__(self, curveDate,
                 timesOrDates,
                 zeroRates,
                 frequencyType=FinFrequencyTypes.ANNUAL,
                 dayCountType=FinDayCountTypes.ACT_ACT_ISDA,
                 interpMethod=FinInterpMethods.FLAT_FORWARDS):
        ''' Create the discount curve from a vector of times and discount
        factors. First date is the curve anchor and first rates should be zero
        as it starts and ends on that date and so has no impact. '''

        # Validate curve
        if len(timesOrDates) == 0:
            raise FinError("Times has zero length")

        if len(timesOrDates) != len(zeroRates):
            raise FinError("Times and Values are not the same")

        if frequencyType not in FinFrequencyTypes:
            raise FinError("Unknown Frequency type " + str(frequencyType))

        if dayCountType not in FinDayCountTypes:
            raise FinError("Unknown Cap Floor DayCountRule type " +
                           str(dayCountType))

        times = []
        if isinstance(timesOrDates[0], float):

            for i in range(0, len(timesOrDates)):
                t = timesOrDates[i]
                if t < 0.0:
                    raise FinError("Times must be > 0.")
                times.append(t)
        elif isinstance(timesOrDates[0], FinDate):

            dc = FinDayCount(dayCountType)
            for i in range(0, len(timesOrDates)):
                t = dc.yearFrac(curveDate, timesOrDates[i])
                if t < 0.0:
                    raise FinError("Times must be > 0.")
                times.append(t)
        else:
            raise FinError("Input timeOrDates must be list of times or dates")

        times = np.array(times)

        if testMonotonicity(times) is False:
            raise FinError("Times are not sorted in increasing order")

        freq = FinFrequency(frequencyType)

        values = []
        for i in range(0, len(timesOrDates)):
            t = times[i]
            r = zeroRates[i]
            df = 1.0 / np.power(1.0 + r/freq, freq * t)
            values.append(df)

        self._curveDate = curveDate
        self._times = times
        self._values = np.array(values)
        self._dayCountType = dayCountType
        self._frequencyType = frequencyType
        self._interpMethod = interpMethod

###############################################################################

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

    def bump(self, bumpSize):
        ''' Calculate the continuous forward rate at the forward date. '''

        times = self._times.copy()
        values = self._values.copy()

        n = len(self._times)
        for i in range(0, n):
            t = times[i]
            values[i] = values[i] * np.exp(-bumpSize*t)

        discCurve = FinDiscountCurve(self._curveDate, times, values,
                                     self._interpMethod)

        return discCurve

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

    def __repr__(self):

        numPoints = len(self._times)
        s = labelToString("TIMES", "DISCOUNT FACTORS")
        for i in range(0, numPoints):
            s += labelToString(self._times[i], self._values[i])

        return s

#######################################################################
