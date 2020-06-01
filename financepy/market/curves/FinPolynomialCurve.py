##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

##########################################################################

from ...market.curves.FinCurve import inputTime
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinHelperFunctions import labelToString

##########################################################################


class FinPolynomialCurve():
    ''' Curve with zero rate of specified frequency parametrised 
    as a cubic polynomial. '''

    def __init__(self, curveDate, coefficients,
                 compoundingType=-1):
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
        r = self.zeroRate(t)
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

##############################################################################

    def __repr__(self):
        ''' Display internal parameters of curve. '''
        s = labelToString("POWER", "COEFFICIENT")
        for i in range(0, len(self._coefficients)):
            s += labelToString(str(i), self._coefficients[i])

        return s

##############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

##############################################################################
