##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import G_DAYS_IN_YEARS

########################################################################################
# TODO: Market calibration (fitting)
########################################################################################


class IborCapVolCurveFn:
    """Class to manage a term structure of caplet volatilities using the
    parametric form suggested by Rebonato (1999)."""

    def __init__(self, curve_dt, a, b, c, d):

        self._curve_dt = curve_dt
        self._a = a
        self._b = b
        self._c = c
        self._d = d

    ###########################################################################

    def cap_floorlet_vol(self, dt):
        """Return the caplet volatility."""

        if isinstance(dt, Date):
            t = (dt - self._curve_dt) / G_DAYS_IN_YEARS
        else:
            t = dt

        vol = (self._a + self._b * t) * np.exp(-self._c * t) + self._d

        if vol < 0.0:
            raise FinError("Negative volatility. Not permitted.")

        return vol


########################################################################################
