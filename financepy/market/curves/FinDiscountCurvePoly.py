##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

###############################################################################

from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinHelperFunctions import labelToString, inputTime
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################


class FinDiscountCurvePoly(FinDiscountCurve):
    ''' Curve with zero rate of specified frequency parametrised
    as a cubic polynomial. '''

    def __init__(self,
                 curveDate,
                 coefficients,
                 compoundingType=-1):
        ''' Create cubic curve from coefficients '''
        self._curveDate = curveDate
        self._coefficients = coefficients
        self._power = len(coefficients) - 1

###############################################################################

    def zeroRate(self, dt):
        ''' Zero rate from polynomial zero curve. '''
        t = inputTime(dt, self)
        zeroRate = 0.0
        for n in range(0, len(self._coefficients)):
            zeroRate += self._coefficients[n] * (t**n)
        return zeroRate

###############################################################################

    def df(self, dt):
        ''' Discount factor from polynomial zero curve.'''
        t = inputTime(dt, self)
        r = self.zeroRate(t)
        return np.exp(-r * t)

###############################################################################

    def fwd(self, dt):
        ''' Continuously compounded forward rate. '''
        t = inputTime(dt, self)
        dzerodt = 0.0
        for n in range(0, len(self._coefficients)):
            dzerodt += n * self._coefficients[n] * (t**(n-1))

        zeroRate = self.zeroRate(t)
        fwdRate = (zeroRate + t*dzerodt)
        return fwdRate

###############################################################################

    def __repr__(self):
        ''' Display internal parameters of curve. '''
        s = labelToString("POWER", "COEFFICIENT")
        for i in range(0, len(self._coefficients)):
            s += labelToString(str(i), self._coefficients[i])

        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
