##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear

##########################################################################
# TODO: Market calibration (fitting)
##########################################################################


class FinIborCapVolCurveFn():
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

    def capFloorletVol(self, dt):
        ''' Return the caplet volatility. '''

        if isinstance(dt, FinDate):
            t = (dt - self._curveDate) / gDaysInYear
            vol = (self._a + self._b*t) * np.exp(-self._c*t) + self._d

        if vol < 0.0:
            raise FinError("Negative volatility. Not permitted.")

        return vol

###############################################################################
