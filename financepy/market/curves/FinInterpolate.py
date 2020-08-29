##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, fabs
from numba import njit, float64, int64
import numpy as np
from ...finutils.FinError import FinError

###############################################################################

from enum import Enum


class FinInterpTypes(Enum):
    LINEAR_ZERO_RATES = 1
    FLAT_FORWARDS = 2
    LINEAR_FORWARDS = 3
    LINEAR_SWAP_RATES = 4

###############################################################################


def interpolate(t: (float, np.ndarray),  # time or array of times
                times: np.ndarray,  # Vector of times on grid
                dfs: np.ndarray,  # Vector of discount factors
                method: int):  # Interpolation method
    ''' Fast interpolation of discount factors at time x given discount factors
    at times provided using one of the methods in the enum FinInterpTypes. The
    value of x can be an array so that the function is vectorised. '''

    if type(t) is float or type(t) is np.float64:
        u = _uinterpolate(t, times, dfs, method)
        return u
    elif type(t) is np.ndarray:
        v = _vinterpolate(t, times, dfs, method)
        return v
    else:
        raise FinError("Unknown input type" + type(t))

###############################################################################


@njit(float64(float64, float64[:], float64[:], int64),
      fastmath=True, cache=True, nogil=True)
def _uinterpolate(t, times, dfs, method):
    ''' Return the interpolated value of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate. '''

    if method == FinInterpTypes.LINEAR_SWAP_RATES.value:
        method = FinInterpTypes.FLAT_FORWARDS.value

    small = 1e-10
    numPoints = times.size

    if t == times[0]:
        return dfs[0]

    i = 0
    while times[i] < t and i < numPoints - 1:
        i = i + 1

    if t > times[i]:
        i = numPoints

    yvalue = 0.0

    ###########################################################################
    # linear interpolation of y(x)
    ###########################################################################

    if method == FinInterpTypes.LINEAR_ZERO_RATES.value:

        if i == 1:
            r1 = -np.log(dfs[i])/times[i]
            r2 = -np.log(dfs[i])/times[i]
            dt = times[i] - times[i-1]
            rvalue = ((times[i]-t)*r1 + (t-times[i-1])*r2)/dt
            yvalue = np.exp(-rvalue*t)
        elif i < numPoints:
            r1 = -np.log(dfs[i-1])/times[i-1]
            r2 = -np.log(dfs[i])/times[i]
            dt = times[i] - times[i-1]
            rvalue = ((times[i]-t)*r1 + (t-times[i-1])*r2)/dt
            yvalue = np.exp(-rvalue*t)
        else:
            r1 = -np.log(dfs[i-1])/times[i-1]
            r2 = -np.log(dfs[i-1])/times[i-1]
            dt = times[i-1] - times[i-2]
            rvalue = ((times[i-1]-t)*r1 + (t-times[i-2])*r2)/dt
            yvalue = np.exp(-rvalue*t)

        return yvalue

    ###########################################################################
    # linear interpolation of log(y(x)) which means the linear interpolation of
    # continuously compounded zero rates in the case of discount curves
    # This is also FLAT FORWARDS
    ###########################################################################

    elif method == FinInterpTypes.FLAT_FORWARDS.value:

        if i == 1:
            rt1 = -np.log(dfs[i-1])
            rt2 = -np.log(dfs[i])
            dt = times[i] - times[i-1]
            rtvalue = ((times[i]-t)*rt1 + (t-times[i-1])*rt2)/dt
            yvalue = np.exp(-rtvalue)
        elif i < numPoints:
            rt1 = -np.log(dfs[i - 1])
            rt2 = -np.log(dfs[i])
            dt = times[i] - times[i-1]
            rtvalue = ((times[i]-t)*rt1 + (t-times[i-1])*rt2)/dt
            yvalue = np.exp(-rtvalue)
        else:
            rt1 = -np.log(dfs[i-2])
            rt2 = -np.log(dfs[i-1])
            dt = times[i-1] - times[i-2]
            rtvalue = ((times[i-1]-t)*rt1 + (t-times[i-2])*rt2)/dt
            yvalue = np.exp(-rtvalue)

        return yvalue

    elif method == FinInterpTypes.LINEAR_FORWARDS.value:

        if i == 1:
            y2 = -log(fabs(dfs[i]) + small)
            yvalue = t * y2 / (times[i] + small)
            yvalue = exp(-yvalue)
        elif i < numPoints:
            # If you get a math domain error it is because you need negativ
            fwd1 = -log(dfs[i-1]/dfs[i-2])/(times[i-1]-times[i-2])
            fwd2 = -log(dfs[i]/dfs[i-1])/(times[i]-times[i-1])
            dt = times[i] - times[i-1]
            fwd = ((times[i]-t)*fwd1 + (t-times[i-1])*fwd2)/dt
            yvalue = dfs[i - 1] * np.exp(-fwd * (t - times[i - 1]))
        else:
            fwd = -np.log(dfs[i-1]/dfs[i-2])/(times[i-1]-times[i-2])
            yvalue = dfs[i-1] * np.exp(-fwd * (t - times[i-1]))

        return yvalue

    else:

        raise FinError("Invalid interpolation scheme.")
        return 0.0

###############################################################################

@njit(float64[:](float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True, nogil=True)
def _vinterpolate(xValues,
                  xvector,
                  dfs,
                  method):
    ''' Return the interpolated values of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate. '''

    n = xValues.size
    yvalues = np.empty(n)
    for i in range(0, n):
        yvalues[i] = _uinterpolate(xValues[i], xvector, dfs, method)

    return yvalues

###############################################################################
