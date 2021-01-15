##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import jit, njit, float64, int64

from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...products.fx.FinFXVanillaOption import FinFXVanillaOption
from ...models.FinModelOptionImpliedDbn import optionImpliedDbn
from ...products.fx.FinFXMktConventions import FinFXATMMethod
from ...products.fx.FinFXMktConventions import FinFXDeltaMethod
from ...finutils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...models.FinModelBlackScholes import FinModelBlackScholes

from ...models.FinModelVolatilityFns import FinVolFunctionTypes
from ...models.FinModelVolatilityFns import volFunctionClark
from ...models.FinModelVolatilityFns import volFunctionBloomberg
from ...models.FinModelVolatilityFns import volFunctionSVI
from ...models.FinModelVolatilityFns import volFunctionSSVI
from ...models.FinModelSABR import volFunctionSABR
from ...models.FinModelSABR import volFunctionSABR_BETA_HALF
from ...models.FinModelSABR import volFunctionSABR_BETA_ONE

from ...models.FinModelVolatilityFns import FinVolFunctionTypes

from ...finutils.FinMath import norminvcdf

from ...models.FinModelBlackScholesAnalytical import bsValue
from ...products.fx.FinFXVanillaOption import fastDelta
from ...finutils.FinDistribution import FinDistribution

from ...finutils.FinSolvers1D import newton_secant
from ...finutils.FinSolversNM import nelder_mead
from ...finutils.FinGlobalTypes import FinSolverTypes

###############################################################################
# ISSUES
# sabr does not fit inverted skew curves like eurjpy 
# problem with initial values ? optimiser can drive vol negative
#
# tried adding a function called gap but it screws up the pdf. Need it to 
# be smooth c3. abandoned for moment. Advise use quintic CLARK5 for best fit
#
# examine other functions for vol
#
# find python version of cg minimiser to apply numba to
###############################################################################
 
###############################################################################
# TODO: Speed up search for strike by providing derivative function to go with
#       delta fit.
###############################################################################

# @njit(fastmath=True, cache=True)
# def _g(K, *args):
#     ''' This is the objective function used in the determination of the FX
#     option implied strike which is computed in the class below. '''

#     s = args[0]
#     t = args[1]
#     rd = args[2]
#     rf = args[3]
#     volatility = args[4]
#     deltaMethodValue = args[5]
#     optionTypeValue = args[6]
#     deltaTarget = args[7]

#     deltaOut = fastDelta(s, t, K, rd, rf,
#                          volatility,
#                          deltaMethodValue,
#                          optionTypeValue)

#     objFn = deltaTarget - deltaOut
#     return objFn

###############################################################################

# @njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
# def _interpolateGap(k, strikes, gaps):

#     if k <= strikes[0]:
#         return 0.0

#     if k >= strikes[-1]:
#         return 0.0

#     index = 0
#     for i in range(1, len(strikes)):
#         if k > strikes[i-1]and k <= strikes[i]:
#             index = i
#             break

#     if index == 0:
#         raise FinError("Value not bracketed")

#     k0 = strikes[index-1]
#     k1 = strikes[index]
#     v0 = gaps[index-1]
#     v1 = gaps[index]
#     v = ((k-k0) * v1 + (k1-k) * v0) / (k1-k0)
#     return v

###############################################################################

@njit(fastmath=True, cache=True)
def _obj(params, *args):
    ''' Return a value that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the volTypeValue at a single time slice only.
    '''

    t = args[0]
    f = args[1]
    strikesGrid = args[2]
    index = args[3]
    volatilityGrid = args[4]
    volTypeValue = args[5]

    tot = 0.0

    numStrikes = len(volatilityGrid)
    numExpiryDates = len(volatilityGrid[0])
 
    for i in range(0, numStrikes):

        k = strikesGrid[i][index]

        fittedVol = volFunction(volTypeValue, params, f, k, t)

        mktVol = volatilityGrid[i][index]

        diff = fittedVol - mktVol
            
        tot += diff**2

    return tot

###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################

def _solveToHorizon(t, f,
                    strikesGrid,
                    timeIndex,
                    volatilityGrid,
                    volTypeValue,
                    xinits,
                    finSolverType):

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    tol = 1e-6

    args = (t, f, strikesGrid, timeIndex, volatilityGrid, volTypeValue)

    # Nelmer-Mead (both SciPy & Numba) is quicker, but occasionally fails 
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if finSolverType == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(_obj, np.array(xinits), 
                               bounds=np.array([[], []]).T, 
                               args=args, tol_f=tol,
                               tol_x=tol, max_iter=1000)
        elif finSolverType == FinSolverTypes.NELDER_MEAD:
            opt = minimize(_obj, xinits, args, method="Nelder-Mead", tol=tol)
            xopt = opt.x
        elif finSolverType == FinSolverTypes.CONJUGATE_GRADIENT:
            opt = minimize(_obj, xinits, args, method="CG", tol=tol)
            xopt = opt.x
    except:
        # If convergence fails try again with CG if necessary
        if finSolverType != FinSolverTypes.CONJUGATE_GRADIENT:
            print('Failed to converge, will try CG')
            opt = minimize(_obj, xinits, args, method="CG", tol=tol)
            xopt = opt.x

    print("t: %9.5f alpha:%9.5f beta: %9.5f rho: %9.5f nu: %9.5f"%
          (t, xopt[0], 0.5, xopt[1], xopt[2]))

    params = np.array(xopt)    
    return params

###############################################################################


@njit(float64(int64, float64[:], float64, float64, float64), 
      cache=True, fastmath=True)
def volFunction(volFunctionTypeValue, params, f, k, t):
    ''' Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. '''

    if volFunctionTypeValue == FinVolFunctionTypes.CLARK.value:
        vol = volFunctionClark(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_ONE.value:
        vol = volFunctionSABR_BETA_ONE(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_HALF.value:
        vol = volFunctionSABR_BETA_HALF(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.BBG.value:
        vol = volFunctionBloomberg(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR.value:
        vol = volFunctionSABR(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.CLARK5.value:
        vol = volFunctionClark(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SVI.value:
        vol = volFunctionSVI(params, f, k, t)
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SSVI.value:
        vol = volFunctionSSVI(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


#@njit(cache=True, fastmath=True)
# def _deltaFit(k, *args):
#     ''' This is the objective function used in the determination of the FX
#     Option implied strike which is computed in the class below. I map it into
#     inverse normcdf space to avoid the flat slope of this function at low vol
#     and high K. It speeds up the code as it allows initial values close to
#     the solution to be used. '''

#     volTypeValue = args[0]
#     s = args[1]
#     t = args[2]
#     rd = args[3]
#     rf = args[4]
#     optionTypeValue = args[5]
#     deltaTypeValue = args[6]
#     inverseDeltaTarget = args[7]
#     params = args[8]
#     strikes = args[9]
#     gaps = args[10]

#     f = s * np.exp((rd-rf)*t)
#     v = volFunction(volTypeValue, params, strikes, gaps, f, k, t)
#     deltaOut = fastDelta(s, t, k, rd, rf, v, deltaTypeValue, optionTypeValue)
#     inverseDeltaOut = norminvcdf(np.abs(deltaOut))
#     invObjFn = inverseDeltaTarget - inverseDeltaOut

#     return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


#@njit(float64(float64, float64, float64, float64, int64, int64, float64,
#              int64, float64, float64[:], float64[:], float64[:]), 
#      fastmath=True)
# def _solveForSmileStrike(s, t, rd, rf,
#                             optionTypeValue,
#                             volatilityTypeValue,
#                             deltaTarget,
#                             deltaMethodValue,
#                             initialGuess,
#                             parameters,
#                             strikes,
#                             gaps):
#     ''' Solve for the strike that sets the delta of the option equal to the
#     target value of delta allowing the volatility to be a function of the
#     strike. '''

#     inverseDeltaTarget = norminvcdf(np.abs(deltaTarget))

#     argtuple = (volatilityTypeValue, s, t, rd, rf, 
#                 optionTypeValue, deltaMethodValue, 
#                 inverseDeltaTarget, 
#                 parameters, strikes, gaps)

#     K = newton_secant(_deltaFit, x0=initialGuess, args=argtuple,
#                       tol=1e-8, maxiter=50)

#     return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


# @njit(float64(float64, float64, float64, float64, int64, float64,
#               int64, float64), fastmath=True)
# def solveForStrike(spotFXRate,
#                    tdel, rd, rf,
#                    optionTypeValue,
#                    deltaTarget,
#                    deltaMethodValue,
#                    volatility):
#     ''' This function determines the implied strike of an FX option
#     given a delta and the other option details. It uses a one-dimensional
#     Newton root search algorith to determine the strike that matches an
#     input volatility. '''

#     # =========================================================================
#     # IMPORTANT NOTE:
#     # =========================================================================
#     # For some delta quotation conventions I can solve for K explicitly.
#     # Note that as I am using the function norminvdelta to calculate the
#     # inverse value of delta, this may not, on a round trip using N(x), give
#     # back the value x as it is calculated to a different number of decimal
#     # places. It should however agree to 6-7 decimal places. Which is OK.
#     # =========================================================================

#     if deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA.value:

#         domDF = np.exp(-rd*tdel)
#         forDF = np.exp(-rf*tdel)

#         if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
#             phi = +1.0
#         else:
#             phi = -1.0

#         F0T = spotFXRate * forDF / domDF
#         vsqrtt = volatility * np.sqrt(tdel)
#         arg = deltaTarget*phi/forDF  # CHECK THIS !!!
#         norminvdelta = norminvcdf(arg)
#         K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
#         return K

#     elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA.value:

#         domDF = np.exp(-rd*tdel)
#         forDF = np.exp(-rf*tdel)

#         if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
#             phi = +1.0
#         else:
#             phi = -1.0

#         F0T = spotFXRate * forDF / domDF
#         vsqrtt = volatility * np.sqrt(tdel)
#         arg = deltaTarget*phi
#         norminvdelta = norminvcdf(arg)
#         K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
#         return K

#     elif deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

#         argtuple = (spotFXRate, tdel, rd, rf, volatility,
#                     deltaMethodValue, optionTypeValue, deltaTarget)

#         K = newton_secant(_g, x0=spotFXRate, args=argtuple,
#                           tol=1e-7, maxiter=50)

#         return K

#     elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

#         argtuple = (spotFXRate, tdel, rd, rf, volatility,
#                     deltaMethodValue, optionTypeValue, deltaTarget)

#         K = newton_secant(_g, x0=spotFXRate, args=argtuple,
#                           tol=1e-7, maxiter=50)

#         return K

#     else:

#         raise FinError("Unknown FinFXDeltaMethod")

###############################################################################


class FinSwaptionVolSurface():
    ''' Class to perform a calibration of a chosen parametrised surface to the
    prices of swaptions at different expiry dates and swap tenors. There is a 
    choice of volatility function from cubic in delta to full SABR and SSVI. 
    Check out FinVolFunctionTypes. Visualising the volatility curve is useful. 
    Also, there is no guarantee that the implied pdf will be positive.'''

    def __init__(self,
                 valueDate: FinDate,
                 expiryDates: (list),
                 fwdSwapRates: (list, np.ndarray),
                 strikeGrid: (np.ndarray),
                 volatilityGrid: (np.ndarray),
                 volatilityFunctionType:FinVolFunctionTypes=FinVolFunctionTypes.SABR,
                 finSolverType:FinSolverTypes=FinSolverTypes.NELDER_MEAD):
        ''' Create the FinSwaptionVolSurface object by passing in market vol 
        data for a list of strikes and expiry dates. '''

        checkArgumentTypes(self.__init__, locals())

        self._valueDate = valueDate

        if len(strikeGrid.shape) != 2:
            raise FinError("Strike grid must be a 2D grid of values")

        if len(volatilityGrid.shape) != 2:
            raise FinError("Volatility grid must be a 2D grid of values")
            
        if len(strikeGrid) != len(volatilityGrid):
            raise FinError("Strike grid and volatility grid must have same size")

        if len(strikeGrid[0]) != len(volatilityGrid[0]):
            raise FinError("Strike grid and volatility grid must have same size")

        if len(expiryDates) != len(volatilityGrid[0]):
            raise FinError("Expiry dates not same size as volatility grid")

        self._numExpiryDates = len(volatilityGrid[0])
        self._numStrikes = len(volatilityGrid)

        self._strikeGrid = strikeGrid
        self._volatilityGrid = volatilityGrid

        self._expiryDates = expiryDates
        self._volatilityFunctionType = volatilityFunctionType

        self._fwdSwapRates = fwdSwapRates

        self._buildVolSurface(finSolverType=finSolverType)

###############################################################################

    def volatilityFromStrikeDate(self, K, expiryDate):
        ''' Interpolates the Black-Scholes volatility from the volatility
        surface given call option strike and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility.'''

        texp = (expiryDate - self._valueDate) / gDaysInYear

        volTypeValue = self._volatilityFunctionType.value

        index0 = 0 # lower index in bracket
        index1 = 0 # upper index in bracket

        numCurves = self._numExpiryDates

        if numCurves == 1:

            index0 = 0
            index1 = 0
            
        # If the time is below first time then assume a flat vol
        elif texp <= self._texp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif texp >= self._texp[-1]:

            index0 = len(self._texp) - 1
            index1 = len(self._texp) - 1

        else: # Otherwise we look for bracketing times and interpolate

            for i in range(1, numCurves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._fwdSwapRates[index0]
        fwd1 = self._fwdSwapRates[index1]
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        vol0 = volFunction(volTypeValue, self._parameters[index0],
                               fwd0, K, t0)

        if index1 != index0:

            vol1 = volFunction(volTypeValue, self._parameters[index1],
                               fwd1, K, t1)

        else:
            
            vol1 = vol0

        # In the expiry time dimension, both volatilities are interpolated 
        # at the same strikes but different deltas.
        vart0 = vol0*vol0*t0
        vart1 = vol1*vol1*t1

        if np.abs(t1-t0) > 1e-6:
            vart = ((texp-t0) * vart1 + (t1-texp) * vart0) / (t1 - t0)

            if vart < 0.0:
                raise FinError("Negative variance.")

            volt = np.sqrt(vart/texp)

        else:
            volt = vol1

        return volt

###############################################################################

    # def deltaToStrike(self, callDelta, expiryDate, deltaMethod):
    #     ''' Interpolates the strike at a delta and expiry date. Linear 
    #     interpolation is used in strike.'''

    #     texp = (expiryDate - self._valueDate) / gDaysInYear

    #     volTypeValue = self._volatilityFunctionType.value

    #     s = self._spotFXRate

    #     if deltaMethod is None:
    #         deltaMethodValue = self._deltaMethod.value
    #     else:
    #         deltaMethodValue = deltaMethod.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     numCurves = self._numVolCurves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if numCurves == 1:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is below first time then assume a flat vol
    #     elif texp <= self._texp[0]:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is beyond the last time then extrapolate with a flat vol
    #     elif texp > self._texp[-1]:
 
    #         index0 = len(self._texp) - 1
    #         index1 = len(self._texp) - 1

    #     else: # Otherwise we look for bracketing times and interpolate

    #         for i in range(1, numCurves):

    #             if texp <= self._texp[i] and texp > self._texp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break

    #     #######################################################################
                
    #     t0 = self._texp[index0]
    #     t1 = self._texp[index1]

    #     initialGuess = self._K_ATM[index0]

    #     K0 = _solveForSmileStrike(s, texp, self._rd[index0], self._rf[index0],
    #                               FinOptionTypes.EUROPEAN_CALL.value,
    #                               volTypeValue, callDelta,
    #                               deltaMethodValue,
    #                               initialGuess,
    #                               self._parameters[index0], 
    #                               self._strikes[index0], 
    #                               self._gaps[index0])

    #     if index1 != index0:

    #         K1 = _solveForSmileStrike(s, texp, 
    #                                   self._rd[index1], 
    #                                   self._rf[index1],
    #                                   FinOptionTypes.EUROPEAN_CALL.value,
    #                                   volTypeValue, callDelta,
    #                                   deltaMethodValue,
    #                                   initialGuess,
    #                                   self._parameters[index1], 
    #                                   self._strikes[index1], 
    #                                   self._gaps[index1])
    #     else:

    #         K1 = K0
            
    #     # In the expiry time dimension, both volatilities are interpolated 
    #     # at the same strikes but different deltas.
 
    #     if np.abs(t1-t0) > 1e-6:

    #         K = ((texp-t0) * K1 + (t1-texp) * K1) / (K1 - K0)

    #     else:

    #         K = K1

    #     return K

###############################################################################
        
    # def volatilityFromDeltaDate(self, callDelta, expiryDate, 
    #                             deltaMethod = None):
    #     ''' Interpolates the Black-Scholes volatility from the volatility
    #     surface given a call option delta and expiry date. Linear interpolation
    #     is done in variance space. The smile strikes at bracketed dates are 
    #     determined by determining the strike that reproduces the provided delta
    #     value. This uses the calibration delta convention, but it can be 
    #     overriden by a provided delta convention. The resulting volatilities 
    #     are then determined for each bracketing expiry time and linear 
    #     interpolation is done in variance space and then converted back to a 
    #     lognormal volatility.'''

    #     texp = (expiryDate - self._valueDate) / gDaysInYear

    #     volTypeValue = self._volatilityFunctionType.value

    #     s = self._spotFXRate

    #     if deltaMethod is None:
    #         deltaMethodValue = self._deltaMethod.value
    #     else:
    #         deltaMethodValue = deltaMethod.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     numCurves = self._numVolCurves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if numCurves == 1:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is below first time then assume a flat vol
    #     elif texp <= self._texp[0]:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is beyond the last time then extrapolate with a flat vol
    #     elif texp > self._texp[-1]:
 
    #         index0 = len(self._texp) - 1
    #         index1 = len(self._texp) - 1

    #     else: # Otherwise we look for bracketing times and interpolate

    #         for i in range(1, numCurves):

    #             if texp <= self._texp[i] and texp > self._texp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break
        
    #     fwd0 = self._F0T[index0]
    #     fwd1 = self._F0T[index1]
                
    #     t0 = self._texp[index0]
    #     t1 = self._texp[index1]

    #     initialGuess = self._K_ATM[index0]

    #     K0 = _solveForSmileStrike(s, texp, self._rd[index0], self._rf[index0],
    #                               FinOptionTypes.EUROPEAN_CALL.value,
    #                               volTypeValue, callDelta,
    #                               deltaMethodValue,
    #                               initialGuess,
    #                               self._parameters[index0], 
    #                               self._strikes[index0], 
    #                               self._gaps[index0])

    #     vol0 = volFunction(volTypeValue, self._parameters[index0],
    #                        self._strikes[index0], self._gaps[index0],
    #                        fwd0, K0, t0)

    #     if index1 != index0:

    #         K1 = _solveForSmileStrike(s, texp, 
    #                                   self._rd[index1], 
    #                                   self._rf[index1],
    #                                   FinOptionTypes.EUROPEAN_CALL.value,
    #                                   volTypeValue, callDelta,
    #                                   deltaMethodValue,
    #                                   initialGuess,
    #                                   self._parameters[index1], 
    #                                   self._strikes[index1], 
    #                                   self._gaps[index1])

    #         vol1 = volFunction(volTypeValue, self._parameters[index1],
    #                            self._strikes[index1], self._gaps[index1],
    #                            fwd1, K1, t1)
    #     else:
    #         vol1 = vol0
            
    #     # In the expiry time dimension, both volatilities are interpolated 
    #     # at the same strikes but different deltas.
    #     vart0 = vol0*vol0*t0
    #     vart1 = vol1*vol1*t1

    #     if np.abs(t1-t0) > 1e-6:

    #         vart = ((texp-t0) * vart1 + (t1-texp) * vart0) / (t1 - t0)
    #         kt = ((texp-t0) * K1 + (t1-texp) * K0) / (t1 - t0)

    #         if vart < 0.0:
    #             raise FinError("Failed interpolation due to negative variance.")

    #         volt = np.sqrt(vart/texp)

    #     else:

    #         volt = vol0
    #         kt = K0

    #     return volt, kt

###############################################################################

    def _buildVolSurface(self, finSolverType=FinSolverTypes.NELDER_MEAD):
        ''' Main function to construct the vol surface. '''

        if self._volatilityFunctionType == FinVolFunctionTypes.CLARK:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_ONE:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_HALF:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.BBG:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
            numParameters = 4
        elif self._volatilityFunctionType == FinVolFunctionTypes.CLARK5:
            numParameters = 5
        elif self._volatilityFunctionType == FinVolFunctionTypes.SVI:
            numParameters = 5
        elif self._volatilityFunctionType == FinVolFunctionTypes.SSVI:
            numParameters = 5
        else:
            print(self._volatilityFunctionType)
            raise FinError("Unknown Model Type")

        numExpiryDates = self._numExpiryDates

        self._parameters = np.zeros([numExpiryDates, numParameters])        
        self._texp = np.zeros(numExpiryDates)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        for i in range(0, numExpiryDates):

            expiryDate = self._expiryDates[i]
            texp = (expiryDate - self._valueDate) / gDaysInYear
            self._texp[i] = texp

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        volTypeValue = self._volatilityFunctionType.value

        xinits = []
        xinit = np.zeros(numParameters)
        xinits.append(xinit)

        for i in range(0, numExpiryDates):

            t = self._texp[i]
            f = self._fwdSwapRates[i]

            res = _solveToHorizon(t, f,
                                  self._strikeGrid,
                                  i,
                                  self._volatilityGrid, 
                                  volTypeValue,
                                  xinits[i],
                                  finSolverType)

            self._parameters[i,:] = res

            xinit = res
            xinits.append(xinit)

###############################################################################

    def checkCalibration(self, verbose: bool, tol: float = 1e-6):
        ''' Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals. '''

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valueDate)
            print("STOCK PRICE:", self._stockPrice)
            print("==========================================================")

        K_dummy = 999

        for i in range(0, self._numExpiryDates):

            expiryDate = self._expiryDates[i]
            print("==========================================================")

            for j in range(0, self._numStrikes):
                
                strike = self._strikeGrid[j][i]

                fittedVol = self.volatilityFromStrikeDate(strike, expiryDate)

                mktVol = self._volatilityGrid[j][i]
                
                diff = fittedVol - mktVol
                
                print("%s %12.3f %7.4f %7.4f %7.5f"% 
                      (expiryDate, strike, 
                       fittedVol*100.0, mktVol*100, diff*100))

        print("==========================================================")
        
###############################################################################

    # def impliedDbns(self, lowS, highS, numIntervals):
    #     ''' Calculate the pdf for each tenor horizon. Returns a list of
    #     FinDistribution objects, one for each tenor horizon. '''

    #     dbns = []

    #     for iTenor in range(0, self._numExpiryDates):

    #         f = self._fwdSwapRates[iTenor]
    #         t = self._texp[iTenor]

    #         dS = (highS - lowS)/ numIntervals

    #         disDF = self._discountCurve._df(t)
    #         divDF = self._dividendCurve._df(t)

    #         r = -np.log(disDF) / t
    #         q = -np.log(divDF) / t

    #         Ks = []
    #         vols = []

    #         for iK in range(0, numIntervals):

    #             k = lowS + iK*dS

    #             vol = volFunction(self._volatilityFunctionType.value, 
    #                               self._parameters[iTenor], 
    #                               f, k, t)

    #             Ks.append(k)
    #             vols.append(vol)

    #         Ks = np.array(Ks)
    #         vols = np.array(vols)

    #         density = optionImpliedDbn(self._stockPrice, t, r, q, Ks, vols)

    #         dbn = FinDistribution(Ks, density)
    #         dbns.append(dbn)

    #     return dbns

###############################################################################

    def plotVolCurves(self):
        ''' Generates a plot of each of the vol curves implied by the market 
        and fitted. '''
        
        volTypeVal = self._volatilityFunctionType.value

        for tenorIndex in range(0, self._numExpiryDates):

            lowK = self._strikeGrid[0][tenorIndex] * 0.9
            highK = self._strikeGrid[-1][tenorIndex] * 1.1

            expiryDate = self._expiryDates[tenorIndex]
            plt.figure()

            ks = []
            fittedVols = []

            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            params = self._parameters[tenorIndex]
            t = self._texp[tenorIndex]
            f = self._fwdSwapRates[tenorIndex]
            
            for i in range(0, numIntervals):

                ks.append(K)
                fittedVol = self.volatilityFromStrikeDate(K, expiryDate) * 100.
                fittedVols.append(fittedVol)
                K = K + dK

            labelStr = "FITTED AT " + str(self._expiryDates[tenorIndex])
            plt.plot(ks, fittedVols, label=labelStr)

            labelStr = "MARKET AT " + str(self._expiryDates[tenorIndex])
            mktVols = self._volatilityGrid[:, tenorIndex] * 100.0
            plt.plot(self._strikeGrid[:, tenorIndex], mktVols, 'o', label=labelStr)

            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = str(self._volatilityFunctionType)
            plt.title(title)
            plt.legend()

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUE DATE", self._valueDate)
        s += labelToString("STOCK PRICE", self._stockPrice)
        s += labelToString("ATM METHOD", self._atmMethod)
        s += labelToString("DELTA METHOD", self._deltaMethod)
        s += labelToString("VOL FUNCTION", self._volatilityFunctionType)

        for i in range(0, self._numExpiryDates):

            s += "\n"

            s += labelToString("EXPIRY DATE", self._expiryDates[i])
            s += labelToString("TIME (YRS)", self._texp[i])
            s += labelToString("FWD FX", self._F0T[i])

            s += labelToString("ATM VOLS", self._atmVols[i]*100.0)

            for j in range(0, self._numStrikes):
                
                k = self._strikes[j]
                vol = self._volGrid[i][j]
                print(expiryDate, k, vol)

        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################
