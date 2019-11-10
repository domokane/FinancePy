# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear

##########################################################################
# TODO
##########################################################################

from ...market.curves.FinCurve import FinCurve, FinCompoundingMethods


class FinFlatCurve(FinCurve):
    ''' A simple discount curve based on a single zero rate with its own
    specified compounding method. Hence the curve is assumed to be flat. '''

    def __init__(self,
                 curveDate,
                 rate,
                 compoundingType=FinCompoundingMethods.CONTINUOUS):

        self._parentType = FinCurve
        self._rate = rate

        if not isinstance(curveDate, FinDate):
            raise FinError("CurveDate is not a date " + str(curveDate))

        self._curveDate = curveDate

        if compoundingType == FinCompoundingMethods.CONTINUOUS:
            self._compoundingType = compoundingType
        elif compoundingType == FinCompoundingMethods.ANNUAL:
            self._compoundingType = compoundingType
        elif compoundingType == FinCompoundingMethods.SEMI_ANNUAL:
            self._compoundingType = compoundingType
        elif compoundingType == FinCompoundingMethods.QUARTERLY:
            self._compoundingType = compoundingType
        elif compoundingType == FinCompoundingMethods.MONTHLY:
            self._compoundingType = compoundingType
        elif compoundingType == FinCompoundingMethods.MONEY_MARKET:
            self._compoundingType = compoundingType
        else:
            raise FinError("Unknown compound method " + str(compoundingType))

##########################################################################

    def zero(self, t):
        ''' Return the zero rate which is simply the curve rate. '''
        return self._rate

##########################################################################

    def fwd(self, t):
        ''' Return the fwd rate which is simply the zero rate. '''
        fwdRate = self._rate
        return fwdRate

##########################################################################

    def df(self, time):
        ''' Return the discount factor based on the compounding approach. '''

        if isinstance(time, FinDate):
            t = (time._excelDate - self._curveDate._excelDate) / gDaysInYear
        else:
            t = time

        if self._compoundingType == FinCompoundingMethods.CONTINUOUS:
            df = np.exp(-self._rate * t)
        elif self._compoundingType == FinCompoundingMethods.ANNUAL:
            df = ((1.0 + self._rate)**(-t))
        elif self._compoundingType == FinCompoundingMethods.SEMI_ANNUAL:
            df = ((1.0 + self._rate / 2)**(-t * 2))
        elif self._compoundingType == FinCompoundingMethods.QUARTERLY:
            df = ((1.0 + self._rate / 4)**(-t * 4))
        elif self._compoundingType == FinCompoundingMethods.MONTHLY:
            df = ((1.0 + self._rate / 12)**(-t * 12))
        elif self._compoundingType == FinCompoundingMethods.MONEY_MARKET:
            df = 1.0 / (1.0 + self._rate * t)
        else:
            raise FinError("Unknown compounding method")

        return df

##########################################################################
