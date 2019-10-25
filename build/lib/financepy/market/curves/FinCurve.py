# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinInterpolate import FinInterpMethods

################################################################################
# TO DO
# Write interpolation scheme
################################################################################

from enum import Enum

class FinCompoundingMethods(Enum):
    CONTINUOUS = 1
    ANNUAL = 2
    SEMI_ANNUAL = 3
    QUARTERLY = 4
    MONTHLY = 5
    MONEY_MARKET = 6

################################################################################
    
class FinCurve():

    def __init__(self,curveDate, interpolationMethod, type):

        self._curveDate = curveDate
        self._type = None
        self._times = None
        self._values = None

    def df(self,t, interpolationMethod = FinInterpMethods.FLAT_FORWARDS):

        if t < self._curveDate:
            raise FinError("FinCurve: time before curve date")

        if type(t) is FinDate:
            t = (t - self._curveDate) / gDaysInYear

        return self.interpolate(t)

###############################################################################

