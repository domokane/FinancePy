# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinInterpolate import interpolate
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinError import FinError
import numpy as np

################################################################################

class FinDiscountCurve():

    def __init__(self, curveDate, times, values, method):

        # Validate curve
        if len(times) < 1:
            raise FinError("Times has zero length")

        if len(times) != len(values):
            raise FinError("Times and Values are not the same")

        num_times = len(times)
        
        if times[0] < 0.0:
            raise FinError("First time is negative")

        for i in range(1,num_times):

            if times[i] <= times[i-1]:
                raise FinError("Times are not sorted in increasing order")
            
        self._curveDate = curveDate
        self._times = np.array(times)
        self._values = np.array(values)
        self._method = method

################################################################################

    def df(self, time):

        if type(time) is FinDate:
            t = (time._excelDate - self._curveDate._excelDate) / gDaysInYear
        else:
            t = time

        return interpolate(t, self._times, self._values, self._method.value)

################################################################################
        
    def survivalProbability(self, time):
        
        if type(time) is FinDate:
            t = (time._excelDate - self._curveDate._excelDate) / gDaysInYear
        else: 
            t = time
            
        return interpolate(t, self._times, self._values, self._method.value)
        
################################################################################

    def dump(self):
        
        print("FinDiscountCurve")
        n = len(self._times)
            
        for i in range(0,n):
            print("Time:",self._times[i], "Value:",self._values[i])

################################################################################