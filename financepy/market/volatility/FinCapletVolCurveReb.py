##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinMath import testMonotonicity
from ...finutils.FinHelperFunctions import labelToString
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinDayCount import FinDayCount, FinDayCountTypes

##########################################################################
# TODO: This is unfinished


class FinCapletVolCurveReb():
    ''' Class to manage a term structure of caplet volatilities using the 
    parametric form suggested by Rebonato (1999). '''

    def __init__(self,
                 curveDate,
                 a, 
                 b, 
                 c,
                 d):

        self._curveDate = curveDate
        self._a = a
        self._b = b
        self._c = c
        self._d = d

###############################################################################

    def volatility(self, date):
        ''' Return the caplet volatility. '''

        if isinstance(dt, FinDate): 
            t = (date - self._curveDate) / gDaysInYear

        vol = (a + b*t) * np.exp(-c*t) + d

        if vol < 0.0:
            raise ValueError("Negative volatility. Not permitted.")

        return vol

###############################################################################
