##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np
from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gSmall


from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline

###############################################################################

from enum import Enum


class FinInterpTypes(Enum):
                FLAT_FWD_RATES = 1
                LINEAR_FWD_RATES = 2
                LINEAR_ZERO_RATES = 4
                FINCUBIC_ZERO_RATES = 7
                NATCUBIC_LOG_DISCOUNT = 8
                NATCUBIC_ZERO_RATES = 9
                PCHIP_ZERO_RATES = 10
                PCHIP_LOG_DISCOUNT = 11

# LINEAR_SWAP_RATES = 3

###############################################################################
# TODO: GET RID OF THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###############################################################################

def interpolate(t: (float, np.ndarray),  # time or array of times
                times: np.ndarray,  # Vector of times on grid
                dfs: np.ndarray,  # Vector of discount factors
                method: int):  # Interpolation method which is value of enum
    ''' Fast interpolation of discount factors at time x given discount factors
    at times provided using one of the methods in the enum FinInterpTypes. The
    value of x can be an array so that the function is vectorised. '''

    if type(t) is float or type(t) is np.float64:
        
        if t < 0.0:
            print(t)
            raise FinError("Interpolate times must all be >= 0")

        u = _uinterpolate(t, times, dfs, method)
        return u
    elif type(t) is np.ndarray:

        if np.any(t < 0.0):
            print(t)
            raise FinError("Interpolate times must all be >= 0")

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

    elif method == FinInterpTypes.FLAT_FWD_RATES.value:

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

    elif method == FinInterpTypes.LINEAR_FWD_RATES.value:

        if i == 1:
            y2 = -np.log(dfs[i] + small)
            yvalue = t * y2 / (times[i] + small)
            yvalue = np.exp(-yvalue)
        elif i < numPoints:
            # If you get a math domain error it is because you need negativ
            fwd1 = -np.log(dfs[i-1]/dfs[i-2])/(times[i-1]-times[i-2])
            fwd2 = -np.log(dfs[i]/dfs[i-1])/(times[i]-times[i-1])
            dt = times[i] - times[i-1]
            fwd = ((times[i]-t)*fwd1 + (t-times[i-1])*fwd2)/dt
            yvalue = dfs[i - 1] * np.exp(-fwd * (t - times[i - 1]))
        else:
            fwd = -np.log(dfs[i-1]/dfs[i-2])/(times[i-1]-times[i-2])
            yvalue = dfs[i-1] * np.exp(-fwd * (t - times[i-1]))

        return yvalue

    else:
        print(method)
        raise FinError("Invalid interpolation scheme.")

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


class FinInterpolator():

    def __init__(self,
                 interpolatorType: FinInterpTypes):
        
        self._interpType = interpolatorType
        self._interpFn = None
        self._times = None
        self._dfs = None
        self._refitCurve = False
        
    ###########################################################################
    
    def fit(self,
            times: np.ndarray,
            dfs: np.ndarray):

        self._times = times
        self._dfs = dfs

        if len(times) == 1:
            return

        if self._interpType == FinInterpTypes.PCHIP_LOG_DISCOUNT:

             logDfs = np.log(self._dfs)
             self._interpFn = PchipInterpolator(self._times, logDfs)

        elif self._interpType  == FinInterpTypes.PCHIP_ZERO_RATES:

             gSmallVector = np.ones(len(self._times)) * gSmall
             zeroRates = -np.log(self._dfs)/(self._times + gSmallVector)

             if self._times[0] == 0.0:
                 zeroRates[0] = zeroRates[1]

             self._interpFn = PchipInterpolator(self._times, zeroRates)

        # if self._interpType == FinInterpTypes.FINCUBIC_LOG_DISCOUNT:

        #     ''' Second derivatives at left is zero and first derivative at
        #     right is clamped to zero. '''
        #     logDfs = np.log(self._dfs)
        #     self._interpFn = CubicSpline(self._times, logDfs,
        #                                  bc_type=((2, 0.0), (1, 0.0)))

        elif self._interpType == FinInterpTypes.FINCUBIC_ZERO_RATES:

            ''' Second derivatives at left is zero and first derivative at
            right is clamped to zero. '''
            gSmallVector = np.ones(len(self._times)) * gSmall
            zeroRates = -np.log(self._dfs)/(self._times + gSmallVector)

            if self._times[0] == 0.0:
                zeroRates[0] = zeroRates[1]

            self._interpFn = CubicSpline(self._times, zeroRates,
                                         bc_type=((2, 0.0), (1, 0.0)))

        elif self._interpType == FinInterpTypes.NATCUBIC_LOG_DISCOUNT:

            ''' Second derivatives are clamped to zero at end points '''
            logDfs = np.log(self._dfs)
            self._interpFn = CubicSpline(self._times, logDfs,
                                         bc_type = 'natural')
    
        elif self._interpType == FinInterpTypes.NATCUBIC_ZERO_RATES:

            ''' Second derivatives are clamped to zero at end points '''
            gSmallVector = np.ones(len(self._times)) * gSmall
            zeroRates = -np.log(self._dfs)/(self._times + gSmallVector)

            if self._times[0] == 0.0:
                zeroRates[0] = zeroRates[1]

            self._interpFn = CubicSpline(self._times, zeroRates,
                                         bc_type = 'natural')

#        elif self._interpType  == FinInterpTypes.LINEAR_LOG_DISCOUNT:
#
#            logDfs = np.log(self._dfs)
#            self._interpFn = interp1d(self._times, logDfs, 
#                                      fill_value="extrapolate")

            
    ###########################################################################

    def interpolate(self,
                    t: float):
        ''' Interpolation of discount factors at time x given discount factors
        at times provided using one of the methods in the enum FinInterpTypes.
        The value of x can be an array so that the function is vectorised. '''

        if self._dfs is None:
            raise FinError("Dfs have not been set.")

        if type(t) is float or type(t) is np.float64:

            if t < 0.0:
                print(t)
                raise FinError("Interpolate times must all be >= 0")

            if np.abs(t) < gSmall:
                return 1.0

            tvec = np.array([t])

        elif type(t) is np.ndarray:

            if np.any(t < 0.0):
                print(t)
                raise FinError("Interpolate times must all be >= 0")

            tvec = t

        else:
            raise FinError("t is not a recognized type")
        
        if self._interpType == FinInterpTypes.PCHIP_LOG_DISCOUNT:
        
            out = np.exp(self._interpFn(tvec))

        elif self._interpType == FinInterpTypes.PCHIP_ZERO_RATES:
    
            out = np.exp(-tvec * self._interpFn(tvec))

        # if self._interpType == FinInterpTypes.FINCUBIC_LOG_DISCOUNT:
        
        #     out = np.exp(self._interpFn(tvec))

        elif self._interpType == FinInterpTypes.FINCUBIC_ZERO_RATES:
        
            out = np.exp(-tvec * self._interpFn(tvec))

        elif self._interpType == FinInterpTypes.NATCUBIC_LOG_DISCOUNT:
        
            out = np.exp(self._interpFn(tvec))

        elif self._interpType == FinInterpTypes.NATCUBIC_ZERO_RATES:
        
            out = np.exp(-tvec * self._interpFn(tvec))
    
#        elif self._interpType == FinInterpTypes.LINEAR_LOG_DISCOUNT:
#    
#            out = np.exp(self._interpFn(tvec))

        else:

            out = _vinterpolate(tvec, self._times, self._dfs,
                                self._interpType.value)

        if type(t) is float or type(t) is np.float64:
            return out[0]
        else:
            return out

###############################################################################
