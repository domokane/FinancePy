# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:52:16 2019

@author: Dominic O'Kane
"""

import numpy as np
from numba import jit

from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear

###############################################################################


def inputFrequency(f):
    if f in [-1, 0, 1, 2, 3, 4, 6, 12]:
        return f
    else:
        raise FinError("Unknown frequency" + str(f))

###############################################################################


#@jit
def inputTime(dt, curve):

    small = 1e-8

    def check(t):
        if t < 0.0:
            raise FinError("Date is before curve value date.")
        elif t < small:
            t = small
        return t

    if isinstance(dt, float):
        t = dt
        return check(t)
    elif isinstance(dt, FinDate):
        t = (dt - curve._curveDate) / gDaysInYear
        return check(t)
    elif isinstance(dt, np.ndarray):
        t = dt
        if np.any(t) < 0:
            raise FinError("Date is before curve value date.")
        t = np.maximum(small, t)
        return t
    else:
        raise FinError("Unknown type.")

###############################################################################
