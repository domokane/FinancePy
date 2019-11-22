# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np

##########################################################################

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinDayCount import FinDayCount
from ...finutils.FinHelperFunctions import inputTime, inputFrequency
from ...market.curves.FinDiscountCurve import FinDiscountCurve

##########################################################################


class FinFlatCurve():
    ''' A trivally simple curve based on a single zero rate with its own
    specified compounding method. Hence the curve is assumed to be flat. '''

##########################################################################

    def __init__(self, curveDate, rate, compoundingFreq=-1):
        ''' Create a FinFlatCurve which requires a curve date. '''
        if not isinstance(curveDate, FinDate):
            raise FinError("CurveDate is not a date " + str(curveDate))

        self._curveDate = curveDate
        self._rate = rate
        self._cmpdFreq = inputFrequency(compoundingFreq)

##########################################################################

    def zeroRate(self, dt, compoundingFreq):
        ''' Return the zero rate which is simply the curve rate. '''
        t = inputTime(dt, self)
        f = inputFrequency(compoundingFreq)

        if f == self._cmpdFreq:
            return self._rate

        if self._cmpdFreq == 0:
            df = 1.0/(1.0 + t*self._rate)
        elif self._cmpdFreq == -1:
            df = np.exp(-self._rate * t)
        else:
            df = ((1.0+self._rate/self._cmpdFreq)**(-t*self._cmpdFreq))

        if f == 0:  # Simple interest
            r = (1.0/df-1.0)/t
        elif f == -1:  # Continuous
            r = -np.log(df) / t
        else:
            r = (df**(-1.0/t) - 1) * f
        return r

##########################################################################

    def fwd(self, dt):
        ''' Return the fwd rate which is simply the zero rate. '''
        fwdRate = self._rate
        return fwdRate

##########################################################################

    def df(self, dt):
        ''' Return the discount factor based on the compounding approach. '''
        t = inputTime(dt, self)
        if self._cmpdFreq == 0:
            df = 1.0/(1.0 + t*self._rate)
        elif self._cmpdFreq == -1:
            df = np.exp(-self._rate * t)
        else:
            df = ((1.0 + self._rate/self._cmpdFreq)**(-t*self._cmpdFreq))
        return df

##########################################################################

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

##########################################################################
