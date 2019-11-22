# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np

##########################################################################

from ...market.curves.FinCurve import FinCompoundingMethods, inputTime
from ...finutils.FinDayCount import FinDayCount

##########################################################################


class FinPolynomialCurve():
    ''' Curve with zero rate of specified frequency parametrised 
    as a cubic polynomial. '''

    def __init__(self, curveDate, coefficients,
                 compoundingType=FinCompoundingMethods.CONTINUOUS):
        ''' Create cubic curve from coefficients '''
        self._curveDate = curveDate
        self._coefficients = coefficients
        self._power = len(coefficients) - 1

    def zeroRate(self, dt):
        ''' Zero rate from polynomial zero curve. '''
        t = inputTime(dt, self)
        zeroRate = 0.0
        for n in range(0, len(self._coefficients)):
            zeroRate += self._coefficients[n] * (t**n)
        return zeroRate

    def df(self, dt):
        ''' Discount factor from polynomial zero curve.'''
        t = inputTime(dt, self)
        r = self.zero(t)
        return np.exp(-r * t)

    def fwd(self, dt):
        ''' Continuously compounded forward rate. '''
        t = inputTime(dt, self)
        dzerodt = 0.0
        for n in range(0, len(self._coefficients)):
            dzerodt += n * self._coefficients[n] * (t**(n-1))

        zeroRate = self.zeroRate(t)
        fwdRate = (zeroRate + t*dzerodt)
        return fwdRate

    def fwdRate(self, date1, date2, dayCountType):
        ''' Calculate the forward rate according to the specified
        day count convention. '''

        if date1 < self._curveDate:
            raise ValueError("Date1 before curve date.")

        if date2 < date1:
            raise ValueError("Date2 must not be before Date1")

        dayCount = FinDayCount(dayCountType)
        yearFrac = dayCount.yearFrac(date1, date2)
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / yearFrac
        return fwd

    def print(self):
        for i in range(0, len(self._coefficients)):
            print("Power %d Coefficient %10.7f" % (i, self._coefficients[i]))

##########################################################################
