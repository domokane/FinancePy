# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""
from math import exp, log, fabs
from numba import njit, float64, int64, jit
import numpy as np

###############################################################################

from enum import Enum


class FinInterpMethods(Enum):
    PIECEWISE_LINEAR = 1
    PIECEWISE_LOG_LINEAR = 2
    FLAT_FORWARDS = 3

###############################################################################


def interpolate(xValue,
                xvector,
                yvector,
                method):

    if type(xValue) is float or type(xValue) is np.float64:
        u = uinterpolate(xValue, xvector, yvector, method)
        return u
    elif type(xValue) is np.ndarray:
        v = vinterpolate(xValue, xvector, yvector, method)
        return v
    else:
        raise ValueError("Unknown input type")

###############################################################################

@njit(float64(float64, float64[:], float64[:], int64),
      fastmath=True, cache=True, nogil=True)
def uinterpolate(xValue,
                 xvector,
                 yvector,
                 method):
    ''' Return the interpolated value of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate. '''

    small = 1e-10
    numPoints = xvector.size

    if xValue == xvector[0]:
        return yvector[0]

    i = 0
    while xvector[i] < xValue and i < numPoints - 1:
        i = i + 1

    if xValue > xvector[i]:
        i = numPoints

    yvalue = 0.0

    ###########################################################################
    # linear interpolation of y(x)
    ###########################################################################

    if method == FinInterpMethods.PIECEWISE_LINEAR.value:

        if i == 0:
            yvalue = yvector[0]
        elif i < numPoints:
            y1 = yvector[i - 1]
            y2 = yvector[i]
            yvalue = (xvector[i] - xValue) * y1 + \
                (xValue - xvector[i - 1]) * y2
            yvalue = yvalue / (xvector[i] - xvector[i - 1])
        else:
            y1 = yvector[i - 2]
            y2 = yvector[i - 1]
            slope = (y2 - y1) / (xvector[i - 1] - xvector[i - 2])
            yvalue = yvector[i - 1] + slope * (xValue - xvector[i - 1])

        return yvalue

    ###########################################################################
    # linear interpolation of log(y(x)) which means the linear interpolation of
    # continuously compounded zero rates in the case of discount curves.
    ###########################################################################

    elif method == FinInterpMethods.PIECEWISE_LOG_LINEAR.value:

        if i == 0:
            y2 = -log(fabs(yvector[i]) + small)
            yvalue = xValue * y2 / (xvector[i] + small)
            yvalue = exp(-yvalue)
        elif i < numPoints:
            y1 = -log(yvector[i - 1] + small)
            y2 = -log(yvector[i] + small)
            yvalue = (xvector[i] - xValue) * y1 + \
                (xValue - xvector[i - 1]) * y2
            yvalue = yvalue / (xvector[i] - xvector[i - 1])
            yvalue = exp(-yvalue)
        else:
            y1 = yvector[i - 2]
            y2 = yvector[i - 1]
            yvalue = -log(y2 / y1) / (xvector[i - 1] - xvector[i - 2])
            yvalue = y2 * exp(-(xValue - xvector[i - 1]) * yvalue)

        return yvalue

    elif method == FinInterpMethods.FLAT_FORWARDS.value:

        if i == 0:
            y2 = -log(fabs(yvector[i]) + small)
            yvalue = xValue * y2 / (xvector[i] + small)
            yvalue = exp(-yvalue)
        elif i < numPoints:
            # If you get a math domain error it is because you need negativ
            fwd = log(yvector[i - 1] / (yvector[i] + small))
            fwd = fwd / (xvector[i] - xvector[i - 1])
            yvalue = yvector[i - 1] * exp(-fwd * (xValue - xvector[i - 1]))
        else:
            y1 = yvector[i - 2]
            y2 = yvector[i - 1]
            fwd = -log(y1 / y2)
            fwd = fwd / (xvector[i - 2] - xvector[i - 1])
            yvalue = y2 * exp(-fwd * (xValue - xvector[i - 1]))

        return yvalue

    else:
        return 0.0

###############################################################################

@njit(float64[:](float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True, nogil=True)
def vinterpolate(xValues,
                 xvector,
                 yvector,
                 method):
    ''' Return the interpolated values of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate. '''

    n = xValues.size
    yvalues = np.empty(n)
    for i in range(0, n):
        yvalues[i] = uinterpolate(xValues[i], xvector, yvector, method)

    return yvalues

###############################################################################
