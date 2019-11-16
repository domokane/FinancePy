# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from math import log
import numpy as np

from ...finutils.FinInterpolate import uinterpolate, interpolate
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount

##########################################################################


class FinDiscountCurve():

    def __init__(self,
                 curveDate,
                 times,
                 values,
                 interpMethod):

        # Validate curve
        if len(times) < 1:
            raise FinError("Times has zero length")

        if len(times) != len(values):
            raise FinError("Times and Values are not the same")

        num_times = len(times)

        if times[0] < 0.0:
            raise FinError("First time is negative")

        for i in range(1, num_times):

            if times[i] <= times[i - 1]:
                raise FinError("Times are not sorted in increasing order")

        self._curveDate = curveDate
        self._times = np.array(times)
        self._values = np.array(values)
        self._interpMethod = interpMethod

##########################################################################

    def df(self, time):

        if type(time) is FinDate:
            t = (time._excelDate - self._curveDate._excelDate) / gDaysInYear
        else:
            t = time

        if np.any(t) < 0:
            raise FinError("Date is before curve value date.")

        z = interpolate(t, self._times, self._values,
                        self._interpMethod.value)

        return z

##########################################################################

    def survivalProbability(self, maturityDate):

        if maturityDate < self._curveDate:
            raise FinError("Date is before curve value date.")

        if type(maturityDate) is FinDate:
            t = (maturityDate - self._curveDate) / gDaysInYear
        else:
            t = maturityDate

        q = interpolate(t, self._times, self._values, self._interpMethod.value)
        return q

###############################################################################

    def zeroRate(self, maturityDate, frequency=None):
        ''' Calculate the continuous compounded zero rate to maturity date. '''

        if type(maturityDate) is FinDate:
            t = (maturityDate - self._curveDate) / gDaysInYear
        else:
            t = maturityDate

        if np.any(t) < 0:
            raise FinError("Date is before curve value date.")

        t = np.maximum(t, 1e-10)
        df = self.df(t)

        if frequency is None:
            zeroRate = -np.log(df) / t
        else:
            zeroRate = (df**(-1.0/t) - 1) * frequency

        return zeroRate

##########################################################################

    def fwdContinuous(self, forwardDate):
        ''' Calculate the continuous forward rate at the forward date. '''

        if forwardDate < self._curveDate:
            raise FinError("Forward Date before curve value date.")

        tau = (forwardDate - self._curveDate) / gDaysInYear
        dt = 0.000001
        df1 = self.df(tau)
        df2 = self.df(tau + dt)
        fwd = log(df1 / df2) / dt
        return fwd

##########################################################################

    def fwdLibor(self, date1, date2, dayCountType):
        ''' Calculate the Libor forward rate according to the corresponding
        day count convention. '''

        if date1 < self._curveDate:
            raise FinError("Date " +
                           str(date1) +
                           " before curve value date " +
                           str(self._curveDate))

        if date2 < date1:
            raise FinError("Date2 before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1, date2)
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / yearFrac
        return fwd

#######################################################################

    def print(self):
        ''' Print the details '''

        numPoints = len(self._times)

        print("TIMES", "DISCOUNT FACTORS")

        for i in range(0, numPoints):
            print("%10.7f" % self._times[i], "%10.7f" % self._values[i])
