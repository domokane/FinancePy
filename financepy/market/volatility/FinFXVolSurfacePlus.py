##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

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

from ...models.FinModelVolatilityFns import volFunctionClark
from ...models.FinModelVolatilityFns import volFunctionBloomberg
from ...models.FinModelVolatilityFns import FinVolFunctionTypes
from ...models.FinModelSABR import volFunctionSABR
from ...models.FinModelSABR import volFunctionSABR_BETA_ONE
from ...models.FinModelSABR import volFunctionSABR_BETA_HALF

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

@njit(fastmath=True, cache=True)
def _g(K, *args):
    ''' This is the objective function used in the determination of the FX
    option implied strike which is computed in the class below. '''

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    volatility = args[4]
    deltaMethodValue = args[5]
    optionTypeValue = args[6]
    deltaTarget = args[7]

    deltaOut = fastDelta(s, t, K, rd, rf,
                         volatility,
                         deltaMethodValue,
                         optionTypeValue)

    objFn = deltaTarget - deltaOut
    return objFn

###############################################################################

@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def _interpolateGap(k, strikes, gaps):

    if k <= strikes[0]:
        return 0.0

    if k >= strikes[-1]:
        return 0.0

    index = 0
    for i in range(1, len(strikes)):
        if k > strikes[i-1]and k <= strikes[i]:
            index = i
            break

    if index == 0:
        raise FinError("Value not bracketed")

    k0 = strikes[index-1]
    k1 = strikes[index]
    v0 = gaps[index-1]
    v1 = gaps[index]
    v = ((k-k0) * v1 + (k1-k) * v0) / (k1-k0)
    return v

###############################################################################
# Do not cache this function
@njit(fastmath=True) #, cache=True)
def _obj(params, *args):
    ''' Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the volTypeValue
    '''

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    K_ATM = args[4]
    atmVol = args[5]

    K_25D_C_MS = args[6]
    K_25D_P_MS = args[7]
    V_25D_MS_target = args[8]
    target25DRRVol = args[9]

    K_10D_C_MS = args[10]
    K_10D_P_MS = args[11]
    V_10D_MS_target = args[12]
    target10DRRVol = args[13]

    deltaMethodValue = args[14]
    volTypeValue = args[15]
    alpha = args[16]

    strikesNULL = np.zeros(1)
    gapsNULL = np.zeros(1)

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atmCurveVol = volFunction(volTypeValue, params, strikesNULL, gapsNULL, 
                                  f, K_ATM, t)

    termATM = (atmVol - atmCurveVol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25D strikes
    ###########################################################################

    if target25DRRVol > -999.0:

        sigma_K_25D_C_MS = volFunction(volTypeValue, params, 
                                           strikesNULL, gapsNULL,
                                           f, K_25D_C_MS, t)
    
        V_25D_C_MS = bsValue(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                             FinOptionTypes.EUROPEAN_CALL.value)
    
        sigma_K_25D_P_MS = volFunction(volTypeValue, params, 
                                           strikesNULL, gapsNULL,
                                           f, K_25D_P_MS, t)
    
        V_25D_P_MS = bsValue(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                             FinOptionTypes.EUROPEAN_PUT.value)
    
        V_25D_MS = V_25D_C_MS + V_25D_P_MS
        term25D_1 = (V_25D_MS - V_25D_MS_target)**2

    else:
        
        term25D_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target25DRRVol > -999.0:

        K_25D_C = _solveForSmileStrike(s, t, rd, rf,
                                          FinOptionTypes.EUROPEAN_CALL.value,
                                          volTypeValue, +0.2500,
                                          deltaMethodValue, K_25D_C_MS,
                                          params, strikesNULL, gapsNULL)
    
        sigma_K_25D_C = volFunction(volTypeValue, params, 
                                        strikesNULL, gapsNULL,
                                        f, K_25D_C, t)
    
        K_25D_P = _solveForSmileStrike(s, t, rd, rf,
                                          FinOptionTypes.EUROPEAN_PUT.value,
                                          volTypeValue, -0.2500,
                                          deltaMethodValue, K_25D_P_MS,
                                          params, strikesNULL, gapsNULL)
    
        sigma_K_25D_P = volFunction(volTypeValue, params, 
                                        strikesNULL, gapsNULL,
                                        f, K_25D_P, t)
    
        sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
        term25D_2 = (sigma_25D_RR - target25DRRVol)**2

    else:
        
        term25D_2 = 0.0

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 10D strikes
    ###########################################################################

    if target10DRRVol > -999.0:

        sigma_K_10D_C_MS = volFunction(volTypeValue, params, 
                                           strikesNULL, gapsNULL,
                                           f, K_10D_C_MS, t)
    
        V_10D_C_MS = bsValue(s, t, K_10D_C_MS, rd, rf, sigma_K_10D_C_MS,
                             FinOptionTypes.EUROPEAN_CALL.value)
    
        sigma_K_10D_P_MS = volFunction(volTypeValue, params, 
                                           strikesNULL, gapsNULL,
                                           f, K_10D_P_MS, t)
    
        V_10D_P_MS = bsValue(s, t, K_10D_P_MS, rd, rf, sigma_K_10D_P_MS,
                             FinOptionTypes.EUROPEAN_PUT.value)
    
        V_10D_MS = V_10D_C_MS + V_10D_P_MS
        term10D_1 = (V_10D_MS - V_10D_MS_target)**2

    else:
        
        term10D_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target10DRRVol > -999.0:

        K_10D_C = _solveForSmileStrike(s, t, rd, rf,
                                          FinOptionTypes.EUROPEAN_CALL.value,
                                          volTypeValue, +0.1000,
                                          deltaMethodValue, K_10D_C_MS,
                                          params, strikesNULL, gapsNULL)
    
        sigma_K_10D_C = volFunction(volTypeValue, params, 
                                        strikesNULL, gapsNULL,
                                        f, K_10D_C, t)
    
        K_10D_P = _solveForSmileStrike(s, t, rd, rf,
                                          FinOptionTypes.EUROPEAN_PUT.value,
                                          volTypeValue, -0.1000,
                                          deltaMethodValue, K_10D_P_MS,
                                          params, strikesNULL, gapsNULL)
    
        sigma_K_10D_P = volFunction(volTypeValue, params, 
                                        strikesNULL, gapsNULL,
                                        f, K_10D_P, t)
    
        sigma_10D_RR = (sigma_K_10D_C - sigma_K_10D_P)
        term10D_2 = (sigma_10D_RR - target10DRRVol)**2

    else:
        
        term10D_2 = 0.0

    ###########################################################################
    # Alpha interpolates between fitting only ATM and 25D when alpha = 0.0 and
    # fitting only ATM and 10D when alpha = 1.0. Equal when alpha = 0.50.
    ###########################################################################

    tot = termATM
    tot = tot + (1.0 - alpha) * (term25D_1 + term25D_2)
    tot = tot + alpha * (term10D_1 + term10D_2)
    return tot

###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################

# THIS FUNCTION IS NOT USED CURRENTLY
@njit(fastmath=True) #, cache=True)
def _objGAP(gaps, *args):
    ''' Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the volTypeValue
    '''

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    K_ATM = args[4]
    atmVol = args[5]

    K_25D_C_MS = args[6]
    K_25D_P_MS = args[7]
    V_25D_MS_target = args[8]
    target25DRRVol = args[9]

    K_10D_C_MS = args[10]
    K_10D_P_MS = args[11]
    V_10D_MS_target = args[12]
    target10DRRVol = args[13]

    deltaMethodValue = args[14]
    volTypeValue = args[15]
    params = args[16]

    strikes = [K_10D_P_MS, K_25D_P_MS, K_ATM, K_25D_C_MS, K_10D_C_MS]
    strikes = np.array(strikes)

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atmCurveVol = volFunction(volTypeValue, params, strikes, gaps, 
                              f, K_ATM, t)

    print("atmCurveVol", atmCurveVol)

    termATM = (atmVol - atmCurveVol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25D strikes
    ###########################################################################

    sigma_K_25D_C_MS = volFunction(volTypeValue, params, strikes, gaps, 
                                       f, K_25D_C_MS, t)

    print("sigma_K_25D_C_MS", sigma_K_25D_C_MS)

    V_25D_C_MS = bsValue(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                         FinOptionTypes.EUROPEAN_CALL.value)

    sigma_K_25D_P_MS = volFunction(volTypeValue, params, strikes, gaps, 
                                       f, K_25D_P_MS, t)

    print("sigma_K_25D_P_MS", sigma_K_25D_P_MS)

    V_25D_P_MS = bsValue(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                         FinOptionTypes.EUROPEAN_PUT.value)

    V_25D_MS = V_25D_C_MS + V_25D_P_MS
    term25D_1 = (V_25D_MS - V_25D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    K_25D_C = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, +0.2500,
                                      deltaMethodValue, K_25D_C_MS,
                                      params, strikes, gaps)

    sigma_K_25D_C = volFunction(volTypeValue, params, strikes, gaps,
                                    f, K_25D_C, t)

    print("sigma_K_25D_C", sigma_K_25D_C)

    K_25D_P = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.2500,
                                      deltaMethodValue, K_25D_P_MS,
                                      params, strikes, gaps)

    sigma_K_25D_P = volFunction(volTypeValue, params, strikes, gaps,
                                    f, K_25D_P, t)

    print("sigma_K_25D_P", sigma_K_25D_P)

    sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
    term25D_2 = (sigma_25D_RR - target25DRRVol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 10D strikes
    ###########################################################################

    sigma_K_10D_C_MS = volFunction(volTypeValue, params, strikes, gaps, 
                                       f, K_10D_C_MS, t)

    print("sigma_K_10D_C_MS", sigma_K_10D_C_MS)

    V_10D_C_MS = bsValue(s, t, K_10D_C_MS, rd, rf, sigma_K_10D_C_MS,
                         FinOptionTypes.EUROPEAN_CALL.value)

    sigma_K_10D_P_MS = volFunction(volTypeValue, params, strikes, gaps,
                                       f, K_10D_P_MS, t)

    print("sigma_K_10D_P_MS", sigma_K_10D_P_MS)

    V_10D_P_MS = bsValue(s, t, K_10D_P_MS, rd, rf, sigma_K_10D_P_MS,
                         FinOptionTypes.EUROPEAN_PUT.value)

    V_10D_MS = V_10D_C_MS + V_10D_P_MS
    term10D_1 = (V_10D_MS - V_10D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    K_10D_C = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, +0.1000,
                                      deltaMethodValue, K_10D_C_MS,
                                      params, strikes, gaps)

    sigma_K_10D_C = volFunction(volTypeValue, params, strikes, gaps, 
                                    f, K_10D_C, t)

    print("SIGMA_K_10D_C", sigma_K_10D_C)

    print("INIT K_10D_P_MS", K_10D_P_MS)

    K_10D_P = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.1000,
                                      deltaMethodValue, K_10D_P_MS,
                                      params, strikes, gaps)

    print("K_10D_P", K_10D_P)
    sigma_K_10D_P = volFunction(volTypeValue, params, strikes, gaps,
                                    f, K_10D_P, t)

    print("SIGMA_K_10D_P", sigma_K_10D_P)

    sigma_10D_RR = (sigma_K_10D_C - sigma_K_10D_P)
    term10D_2 = (sigma_10D_RR - target10DRRVol)**2

    ###########################################################################
    # Alpha interpolates between fitting only ATM and 25D when alpha = 0.0 and
    # fitting only ATM and 10D when alpha = 1.0. Equal when alpha = 0.50.
    ###########################################################################

    tot = termATM
    tot = tot + (term25D_1 + term25D_2)
    tot = tot + (term10D_1 + term10D_2)
    return tot

###############################################################################

def _solveToHorizon(s, t, rd, rf,
                       K_ATM, atmVol,
                       ms25DVol, rr25DVol,
                       ms10DVol, rr10DVol,
                       deltaMethodValue, volTypeValue,
                       alpha,
                       xinits,
                       ginits,
                       finSolverType,
                       tol):

    ###########################################################################
    # Determine the price of a market strangle from market strangle
    # Need to price a call and put that agree with market strangle
    ###########################################################################

    use10D = True
    use25D = True

    if ms25DVol == -999.0:
        use25D = False

    if ms10DVol == -999.0:
        use10D = False

    if use25D is True:

        vol_25D_MS = atmVol + ms25DVol
    
        K_25D_C_MS = solveForStrike(s, t, rd, rf,
                                    FinOptionTypes.EUROPEAN_CALL.value,
                                    +0.2500,
                                    deltaMethodValue,
                                    vol_25D_MS)
    
        K_25D_P_MS = solveForStrike(s, t, rd, rf,
                                    FinOptionTypes.EUROPEAN_PUT.value,
                                    -0.2500,
                                    deltaMethodValue,
                                    vol_25D_MS)
    
        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        V_25D_C_MS = bsValue(s, t, K_25D_C_MS, rd, rf, vol_25D_MS,
                             FinOptionTypes.EUROPEAN_CALL.value)
    
        V_25D_P_MS = bsValue(s, t, K_25D_P_MS, rd, rf, vol_25D_MS,
                             FinOptionTypes.EUROPEAN_PUT.value)
    
        # Market price of strangle in the domestic currency
        V_25D_MS = V_25D_C_MS + V_25D_P_MS

    else:
        
        vol_25D_MS = -999.0
        K_25D_C_MS = 0.0
        K_25D_P_MS = 0.0
        V_25D_C_MS = 0.0
        V_25D_P_MS = 0.0
        V_25D_MS = 0.0

    ###########################################################################

    if use10D is True:

        vol_10D_MS = atmVol + ms10DVol
    
        K_10D_C_MS = solveForStrike(s, t, rd, rf,
                                    FinOptionTypes.EUROPEAN_CALL.value,
                                    +0.1000,
                                    deltaMethodValue,
                                    vol_10D_MS)
    
        K_10D_P_MS = solveForStrike(s, t, rd, rf,
                                    FinOptionTypes.EUROPEAN_PUT.value,
                                    -0.1000,
                                    deltaMethodValue,
                                    vol_10D_MS)
    
        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        V_10D_C_MS = bsValue(s, t, K_10D_C_MS, rd, rf, vol_10D_MS,
                             FinOptionTypes.EUROPEAN_CALL.value)
    
        V_10D_P_MS = bsValue(s, t, K_10D_P_MS, rd, rf, vol_10D_MS,
                             FinOptionTypes.EUROPEAN_PUT.value)
    
        # Market price of strangle in the domestic currency
        V_10D_MS = V_10D_C_MS + V_10D_P_MS

    else:

        vol_10D_MS = -999.0
        K_10D_C_MS = 0.0
        K_10D_P_MS = 0.0
        V_10D_C_MS = 0.0
        V_10D_P_MS = 0.0
        V_10D_MS = 0.0

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    # tol = 1e-8

    args = (s, t, rd, rf,
            K_ATM, atmVol,
            K_25D_C_MS, K_25D_P_MS, V_25D_MS, rr25DVol,
            K_10D_C_MS, K_10D_P_MS, V_10D_MS, rr10DVol,
            deltaMethodValue, volTypeValue, alpha)

    # Nelmer-Mead (both SciPy & Numba) is quicker, but occasionally fails 
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if finSolverType == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(_obj, np.array(xinits), 
                               bounds=np.array([[], []]).T, args=args, tol_f=tol,
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

    params = np.array(xopt)

    strikes = [K_10D_P_MS, K_25D_P_MS, K_ATM, K_10D_C_MS, K_25D_C_MS]
    strikes = np.array(strikes)
    gaps = np.zeros(5)
    
    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    if 1==0:

        tol = 1e-12
    
        args = (s, t, rd, rf,
                K_ATM, atmVol,
                K_25D_C_MS, K_25D_P_MS, V_25D_MS, rr25DVol,
                K_10D_C_MS, K_10D_P_MS, V_10D_MS, rr10DVol,
                deltaMethodValue, volTypeValue, params)
    
        opt = minimize(_objGAP, ginits, args, method="Nelder-Mead", tol=tol)
        xopt = opt.x
        gaps = np.array(xopt)
        
        print("SOLVED")

# Removed this as it causes discontinuity
#    f = s * np.exp((rd-rf)*t)
#    interpATMVol = volFunction(volTypeValue, params,
#                                   strikes, gaps, f, K_ATM, t)

#    diff = atmVol - interpATMVol
#    gaps[2] = diff
    
    ###########################################################################

    if use25D is False:
        K_25D_C_MS = K_ATM
        K_25D_P_MS = K_ATM

    K_25D_C = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, +0.2500,
                                      deltaMethodValue, K_25D_C_MS,
                                      params, strikes, gaps)

    K_25D_P = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.2500,
                                      deltaMethodValue, K_25D_P_MS,
                                      params, strikes, gaps)

    if use10D is False:
        K_10D_C_MS = K_ATM
        K_10D_P_MS = K_ATM

    K_10D_C = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, +0.1000,
                                      deltaMethodValue, K_10D_C_MS,
                                      params, strikes, gaps)

    K_10D_P = _solveForSmileStrike(s, t, rd, rf,
                                      FinOptionTypes.EUROPEAN_PUT.value,
                                      volTypeValue, -0.1000,
                                      deltaMethodValue, K_10D_P_MS,
                                      params, strikes, gaps)

    return (params, strikes, gaps,
            K_25D_C_MS, K_25D_P_MS, K_25D_C, K_25D_P,
            K_10D_C_MS, K_10D_P_MS, K_10D_C, K_10D_P)

###############################################################################


@njit(float64(int64, float64[:], float64[:], float64[:], 
              float64, float64, float64), cache=True, fastmath=True)
def volFunction(volFunctionTypeValue, params, strikes, gaps, f, k, t):
    ''' Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. '''

#    print("volFunction", volFunctionTypeValue)

    if len(strikes) == 1:
        gapK = 0.0
    else:
        gapK = _interpolateGap(k, strikes, gaps)

    if volFunctionTypeValue == FinVolFunctionTypes.CLARK.value:
        vol = volFunctionClark(params, f, k, t) + gapK
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR.value:
        vol = volFunctionSABR(params, f, k, t)  + gapK
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_HALF.value:
        vol = volFunctionSABR_BETA_HALF(params, f, k, t)  + gapK
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.SABR_BETA_ONE.value:
        vol = volFunctionSABR_BETA_ONE(params, f, k, t)  + gapK
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.BBG.value:
        vol = volFunctionBloomberg(params, f, k, t)  + gapK
        return vol
    elif volFunctionTypeValue == FinVolFunctionTypes.CLARK5.value:
        vol = volFunctionClark(params, f, k, t)  + gapK
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


@njit(cache=True, fastmath=True)
def _deltaFit(k, *args):
    ''' This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol
    and high K. It speeds up the code as it allows initial values close to
    the solution to be used. '''

    volTypeValue = args[0]
    s = args[1]
    t = args[2]
    rd = args[3]
    rf = args[4]
    optionTypeValue = args[5]
    deltaTypeValue = args[6]
    inverseDeltaTarget = args[7]
    params = args[8]
    strikes = args[9]
    gaps = args[10]

    f = s * np.exp((rd-rf)*t)
    v = volFunction(volTypeValue, params, strikes, gaps, f, k, t)
    deltaOut = fastDelta(s, t, k, rd, rf, v, deltaTypeValue, optionTypeValue)
    inverseDeltaOut = norminvcdf(np.abs(deltaOut))
    invObjFn = inverseDeltaTarget - inverseDeltaOut

#    print(k, f, v, deltaOut, invObjFn)

    return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, int64, float64,
              int64, float64, float64[:], float64[:], float64[:]), 
      fastmath=True)
def _solveForSmileStrike(s, t, rd, rf,
                            optionTypeValue,
                            volatilityTypeValue,
                            deltaTarget,
                            deltaMethodValue,
                            initialGuess,
                            parameters,
                            strikes,
                            gaps):
    ''' Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike. '''

    inverseDeltaTarget = norminvcdf(np.abs(deltaTarget))

    argtuple = (volatilityTypeValue, s, t, rd, rf, 
                optionTypeValue, deltaMethodValue, 
                inverseDeltaTarget, 
                parameters, strikes, gaps)

    K = newton_secant(_deltaFit, x0=initialGuess, args=argtuple,
                      tol=1e-8, maxiter=50)

    return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, float64,
              int64, float64), fastmath=True)
def solveForStrike(spotFXRate,
                   tdel, rd, rf,
                   optionTypeValue,
                   deltaTarget,
                   deltaMethodValue,
                   volatility):
    ''' This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorith to determine the strike that matches an
    input volatility. '''

    # =========================================================================
    # IMPORTANT NOTE:
    # =========================================================================
    # For some delta quotation conventions I can solve for K explicitly.
    # Note that as I am using the function norminvdelta to calculate the
    # inverse value of delta, this may not, on a round trip using N(x), give
    # back the value x as it is calculated to a different number of decimal
    # places. It should however agree to 6-7 decimal places. Which is OK.
    # =========================================================================

    if deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA.value:

        domDF = np.exp(-rd*tdel)
        forDF = np.exp(-rf*tdel)

        if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spotFXRate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = deltaTarget*phi/forDF  # CHECK THIS !!!
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
        return K

    elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA.value:

        domDF = np.exp(-rd*tdel)
        forDF = np.exp(-rf*tdel)

        if optionTypeValue == FinOptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spotFXRate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = deltaTarget*phi
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
        return K

    elif deltaMethodValue == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

        argtuple = (spotFXRate, tdel, rd, rf, volatility,
                    deltaMethodValue, optionTypeValue, deltaTarget)

        K = newton_secant(_g, x0=spotFXRate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    elif deltaMethodValue == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

        argtuple = (spotFXRate, tdel, rd, rf, volatility,
                    deltaMethodValue, optionTypeValue, deltaTarget)

        K = newton_secant(_g, x0=spotFXRate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    else:

        raise FinError("Unknown FinFXDeltaMethod")

###############################################################################


class FinFXVolSurfacePlus():
    ''' Class to perform a calibration of a chosen parametrised surface to the
    prices of FX options at different strikes and expiry tenors. The
    calibration inputs are the ATM and 25 and 10 Delta volatilities in terms of
    the market strangle amd risk reversals. There is a choice of volatility
    function from cubic in delta to full SABR. Check out FinVolFunctionTypes.
    Parameter alpha [0,1] is used to interpolate between fitting only 25D when
    alpha=0 to fitting only 10D when alpha=1.0. Alpha=0.5 assigns equal weights
    A vol function with more parameters will give a better fit. Of course. But 
    it might also overfit. Visualising the volatility curve is useful. Also, 
    there is no guarantee that the implied pdf will be positive.'''

    def __init__(self,
                 valueDate: FinDate,
                 spotFXRate: float,
                 currencyPair: str,
                 notionalCurrency: str,
                 domDiscountCurve: FinDiscountCurve,
                 forDiscountCurve: FinDiscountCurve,
                 tenors: (list),
                 atmVols: (list, np.ndarray),
                 mktStrangle25DeltaVols: (list, np.ndarray),
                 riskReversal25DeltaVols: (list, np.ndarray),
                 mktStrangle10DeltaVols: (list, np.ndarray),
                 riskReversal10DeltaVols: (list, np.ndarray),
                 alpha: float,
                 atmMethod:FinFXATMMethod=FinFXATMMethod.FWD_DELTA_NEUTRAL,
                 deltaMethod:FinFXDeltaMethod=FinFXDeltaMethod.SPOT_DELTA,
                 volatilityFunctionType:FinVolFunctionTypes=FinVolFunctionTypes.CLARK,
                 finSolverType:FinSolverTypes=FinSolverTypes.NELDER_MEAD,
                 tol:float=1e-8):
        ''' Create the FinFXVolSurfacePlus object by passing in market vol data
        for ATM, 25 Delta and 10 Delta strikes. The alpha weight shifts the
        fitting between 25D and 10D. Alpha = 0.0 is 100% 25D while alpha = 1.0
        is 100% 10D. An alpha of 0.50 is equally weighted. '''

        # I want to allow Nones for some of the market inputs
        if mktStrangle10DeltaVols is None:
            mktStrangle10DeltaVols = []

        if riskReversal10DeltaVols is None:
            riskReversal10DeltaVols = []

        if mktStrangle25DeltaVols is None:
            mktStrangle25DeltaVols = []

        if riskReversal25DeltaVols is None:
            riskReversal25DeltaVols = []

        checkArgumentTypes(self.__init__, locals())

        self._valueDate = valueDate
        self._spotFXRate = spotFXRate
        self._currencyPair = currencyPair

        if len(currencyPair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._forName = self._currencyPair[0:3]
        self._domName = self._currencyPair[3:6]

        self._notionalCurrency = notionalCurrency
        self._domDiscountCurve = domDiscountCurve
        self._forDiscountCurve = forDiscountCurve
        self._numVolCurves = len(tenors)
        self._tenors = tenors

        if len(atmVols) != self._numVolCurves:
            raise FinError("Number ATM vols must equal number of tenors")

        self._atmVols = np.array(atmVols)/100.0

        self._useMS25DVol = True
        self._useRR25DVol = True
        self._useMS10DVol = True
        self._useRR10DVol = True

        # Some of these can be missing which is signified by length zero
        n = len(mktStrangle25DeltaVols)

        if n != self._numVolCurves and n!= 0:
            raise FinError("Number MS25D vols must equal number of tenors")

        if n == 0:
            self._useMS25DVol = False

        n = len(riskReversal25DeltaVols)

        if n != self._numVolCurves and n!= 0:
            raise FinError("Number RR25D vols must equal number of tenors")

        if n == 0:
            self._useRR25DVol = False

        n = len(mktStrangle10DeltaVols)

        if n != self._numVolCurves and n!= 0:
            raise FinError("Number MS10D vols must equal number of tenors")

        if n == 0:
            self._useMS10DVol = False

        n = len(riskReversal10DeltaVols)

        if n != self._numVolCurves and n!= 0:
            raise FinError("Number RR10D vols must equal number of tenors")

        if n == 0:
            self._useRR10DVol = False

        if self._useMS10DVol != self._useRR10DVol:
            raise FinError("You must provide both 10D RR + 10D MS or neither")

        if self._useMS25DVol != self._useRR25DVol:
            raise FinError("You must provide both 25D RR + 25D MS or neither")

        if self._useMS10DVol is False and self._useMS25DVol is False:
            raise FinError("No MS and RR. You must provide 10D or 25D MS + RR.")
            
        self._mktStrangle25DeltaVols = np.array(mktStrangle25DeltaVols)/100.0
        self._riskReversal25DeltaVols = np.array(riskReversal25DeltaVols)/100.0
        self._mktStrangle10DeltaVols = np.array(mktStrangle10DeltaVols)/100.0
        self._riskReversal10DeltaVols = np.array(riskReversal10DeltaVols)/100.0

        if alpha < 0.0 or alpha > 1.0:
            raise FinError("Alpha must be between 0.0 and 1.0")

        self._alpha = alpha
        
        self._atmMethod = atmMethod
        self._deltaMethod = deltaMethod

        if self._deltaMethod == FinFXDeltaMethod.SPOT_DELTA:
            self._deltaMethodString = "pips_spot_delta"
        elif self._deltaMethod == FinFXDeltaMethod.FORWARD_DELTA:
            self._deltaMethodString = "pips_fwd_delta"
        elif self._deltaMethod == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ:
            self._deltaMethodString = "pct_spot_delta_prem_adj"
        elif self._deltaMethod == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ:
            self._deltaMethodString = "pct_fwd_delta_prem_adj"
        else:
            raise FinError("Unknown Delta Type")

        self._volatilityFunctionType = volatilityFunctionType
        self._tenorIndex = 0

        self._expiryDates = []
        for i in range(0, self._numVolCurves):
            expiryDate = valueDate.addTenor(tenors[i])
            self._expiryDates.append(expiryDate)

        self._buildVolSurface(finSolverType=finSolverType, tol=tol)

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

        numCurves = self._numVolCurves

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

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        vol0 = volFunction(volTypeValue, self._parameters[index0],
                               self._strikes[index0], self._gaps[index0],
                               fwd0, K, t0)

        if index1 != index0:

            vol1 = volFunction(volTypeValue, self._parameters[index1],
                               self._strikes[index1], self._gaps[index1],
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

    def deltaToStrike(self, callDelta, expiryDate, deltaMethod):
        ''' Interpolates the strike at a delta and expiry date. Linear 
        interpolation is used in strike.'''

        texp = (expiryDate - self._valueDate) / gDaysInYear

        volTypeValue = self._volatilityFunctionType.value

        s = self._spotFXRate

        if deltaMethod is None:
            deltaMethodValue = self._deltaMethod.value
        else:
            deltaMethodValue = deltaMethod.value

        index0 = 0 # lower index in bracket
        index1 = 0 # upper index in bracket

        numCurves = self._numVolCurves

        # If there is only one time horizon then assume flat vol to this time
        if numCurves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif texp <= self._texp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif texp > self._texp[-1]:
 
            index0 = len(self._texp) - 1
            index1 = len(self._texp) - 1

        else: # Otherwise we look for bracketing times and interpolate

            for i in range(1, numCurves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        #######################################################################
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._K_ATM[index0]

        K0 = _solveForSmileStrike(s, texp, self._rd[index0], self._rf[index0],
                                  FinOptionTypes.EUROPEAN_CALL.value,
                                  volTypeValue, callDelta,
                                  deltaMethodValue,
                                  initialGuess,
                                  self._parameters[index0], 
                                  self._strikes[index0], 
                                  self._gaps[index0])

        if index1 != index0:

            K1 = _solveForSmileStrike(s, texp, 
                                      self._rd[index1], 
                                      self._rf[index1],
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, callDelta,
                                      deltaMethodValue,
                                      initialGuess,
                                      self._parameters[index1], 
                                      self._strikes[index1], 
                                      self._gaps[index1])
        else:

            K1 = K0
            
        # In the expiry time dimension, both volatilities are interpolated 
        # at the same strikes but different deltas.
 
        if np.abs(t1-t0) > 1e-6:

            K = ((texp-t0) * K1 + (t1-texp) * K0) / (K1 - K0)

        else:

            K = K1

        return K

###############################################################################
        
    def volatilityFromDeltaDate(self, callDelta, expiryDate, 
                                deltaMethod = None):
        ''' Interpolates the Black-Scholes volatility from the volatility
        surface given a call option delta and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility.'''

        texp = (expiryDate - self._valueDate) / gDaysInYear

        volTypeValue = self._volatilityFunctionType.value

        s = self._spotFXRate

        if deltaMethod is None:
            deltaMethodValue = self._deltaMethod.value
        else:
            deltaMethodValue = deltaMethod.value

        index0 = 0 # lower index in bracket
        index1 = 0 # upper index in bracket

        numCurves = self._numVolCurves

        # If there is only one time horizon then assume flat vol to this time
        if numCurves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif texp <= self._texp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif texp > self._texp[-1]:
 
            index0 = len(self._texp) - 1
            index1 = len(self._texp) - 1

        else: # Otherwise we look for bracketing times and interpolate

            for i in range(1, numCurves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break
        
        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]
                
        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._K_ATM[index0]

        K0 = _solveForSmileStrike(s, texp, self._rd[index0], self._rf[index0],
                                  FinOptionTypes.EUROPEAN_CALL.value,
                                  volTypeValue, callDelta,
                                  deltaMethodValue,
                                  initialGuess,
                                  self._parameters[index0], 
                                  self._strikes[index0], 
                                  self._gaps[index0])

        vol0 = volFunction(volTypeValue, self._parameters[index0],
                           self._strikes[index0], self._gaps[index0],
                           fwd0, K0, t0)

        if index1 != index0:

            K1 = _solveForSmileStrike(s, texp, 
                                      self._rd[index1], 
                                      self._rf[index1],
                                      FinOptionTypes.EUROPEAN_CALL.value,
                                      volTypeValue, callDelta,
                                      deltaMethodValue,
                                      initialGuess,
                                      self._parameters[index1], 
                                      self._strikes[index1], 
                                      self._gaps[index1])

            vol1 = volFunction(volTypeValue, self._parameters[index1],
                               self._strikes[index1], self._gaps[index1],
                               fwd1, K1, t1)
        else:
            vol1 = vol0
            
        # In the expiry time dimension, both volatilities are interpolated 
        # at the same strikes but different deltas.
        vart0 = vol0*vol0*t0
        vart1 = vol1*vol1*t1

        if np.abs(t1-t0) > 1e-6:

            vart = ((texp-t0) * vart1 + (t1-texp) * vart0) / (t1 - t0)
            kt = ((texp-t0) * K1 + (t1-texp) * K0) / (t1 - t0)

            if vart < 0.0:
                raise FinError("Failed interpolation due to negative variance.")

            volt = np.sqrt(vart/texp)

        else:

            volt = vol0
            kt = K0

        return volt, kt

###############################################################################

    def _buildVolSurface(self, finSolverType=FinSolverTypes.NELDER_MEAD, tol=1e-8):
        ''' Main function to construct the vol surface. '''

        s = self._spotFXRate
        numVolCurves = self._numVolCurves

        if self._volatilityFunctionType == FinVolFunctionTypes.CLARK:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
            numParameters = 4
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_ONE:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_HALF:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.BBG:
            numParameters = 3
        elif self._volatilityFunctionType == FinVolFunctionTypes.CLARK5:
            numParameters = 5
        else:
            print(self._volatilityFunctionType)
            raise FinError("Unknown Model Type")

        self._parameters = np.zeros([numVolCurves, numParameters])

        numStrikes = 5
        self._strikes = np.zeros([numVolCurves, numStrikes])
        self._gaps = np.zeros([numVolCurves, numStrikes])
        
        self._texp = np.zeros(numVolCurves)

        self._F0T = np.zeros(numVolCurves)
        self._rd = np.zeros(numVolCurves)
        self._rf = np.zeros(numVolCurves)
        self._K_ATM = np.zeros(numVolCurves)
        self._deltaATM = np.zeros(numVolCurves)

        self._K_25D_C = np.zeros(numVolCurves)
        self._K_25D_P = np.zeros(numVolCurves)
        self._K_25D_C_MS = np.zeros(numVolCurves)
        self._K_25D_P_MS = np.zeros(numVolCurves)
        self._V_25D_MS = np.zeros(numVolCurves)

        self._K_10D_C = np.zeros(numVolCurves)
        self._K_10D_P = np.zeros(numVolCurves)
        self._K_10D_C_MS = np.zeros(numVolCurves)
        self._K_10D_P_MS = np.zeros(numVolCurves)
        self._V_10D_MS = np.zeros(numVolCurves)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################
        spotDate = self._valueDate

        for i in range(0, numVolCurves):

            expiryDate = self._expiryDates[i]
            texp = (expiryDate - spotDate) / gDaysInYear

            domDF = self._domDiscountCurve._df(texp)
            forDF = self._forDiscountCurve._df(texp)
            f = s * forDF/domDF

            self._texp[i] = texp
            self._rd[i] = -np.log(domDF) / texp
            self._rf[i] = -np.log(forDF) / texp
            self._F0T[i] = f

            atmVol = self._atmVols[i]

            # This follows exposition in Clarke Page 52
            if self._atmMethod == FinFXATMMethod.SPOT:
                self._K_ATM[i] = s
            elif self._atmMethod == FinFXATMMethod.FWD:
                self._K_ATM[i] = f
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL:
                self._K_ATM[i] = f * np.exp(atmVol*atmVol*texp/2.0)
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ:
                self._K_ATM[i] = f * np.exp(-atmVol*atmVol*texp/2.0)
            else:
                raise FinError("Unknown Delta Type")

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        ginit = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        xinits = []
        ginits = []
        for i in range(0, numVolCurves):

            atmVol = self._atmVols[i]

            if self._useMS25DVol > 0:
                ms25 = self._mktStrangle25DeltaVols[i]
            else:
                ms25 = 0.0

            if self._useRR25DVol > 0:
                rr25 = self._riskReversal25DeltaVols[i]
            else:
                rr25 = 0.0

            if self._useMS10DVol > 0:
                ms10 = self._mktStrangle10DeltaVols[i]
            else:
                ms10 = 0.0

            if self._useRR10DVol > 0:
                rr10 = self._riskReversal10DeltaVols[i]
            else:
                rr10 = 0.0

            # https://quantpie.co.uk/fx/fx_rr_str.php
            s25 = atmVol + ms25 + rr25/2.0 # 25D Call
            s50 = atmVol                   # ATM
            s75 = atmVol + ms25 - rr25/2.0 # 25D Put (75D Call)

            s10 = atmVol + ms10 + rr10/2.0 # 10D Call
            s50 = atmVol                   # ATM
            s90 = atmVol + ms10 - rr10/2.0 # 10D Put (90D Call)

            if self._volatilityFunctionType == FinVolFunctionTypes.CLARK:

                # Our preference is to fit to the 10D wings first
                if self._useMS10DVol is False:
                    # Fit to 25D
                    c0 = np.log(atmVol)
                    c1 = 2.0 * np.log(s75/s25)
                    c2 = 8.0 * np.log(s25*s75/atmVol/atmVol)
                    xinit = [c0, c1, c2]
                else:
                    # Fit to 10D
                    c0 = np.log(atmVol)
                    c1 = np.log(s90/s10) / 0.80
                    c2 = np.log(s10*s90/atmVol/atmVol) / 0.32
                    xinit = [c0, c1, c2]

            elif self._volatilityFunctionType == FinVolFunctionTypes.SABR:
                # SABR parameters are alpha, nu, rho
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0
                rho = -0.112
                nu = 0.817

                xinit = [alpha, beta, rho, nu]

            elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_ONE:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0 # FIXED
                rho = -0.112
                nu = 0.817

                xinit = [alpha, rho, nu]

            elif self._volatilityFunctionType == FinVolFunctionTypes.SABR_BETA_HALF:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 0.50 # FIXED
                rho = -0.112
                nu = 0.817
 
                xinit = [alpha, rho, nu]

            elif self._volatilityFunctionType == FinVolFunctionTypes.BBG:

                # Our preference is to fit to the 10D wings first
                if self._useMS10DVol is False:
                    # BBG Params if we fit to 25D
                    a = 8.0*s75-16.0*s50+8.0*s25
                    b = -6.0*s75+16.0*s50-10.0*s25
                    c = s75-3.0*s50+3.0*s25
                else:
                    # BBG Params if we fit to 10D
                    a = (25.0*s90 - 50.0*s50 + 25.0*s10) / 8.0
                    b = (-15.0*s90 + 50.0*s50 - 35.0*s10) / 8.0
                    c = (5.0*s90 - 18.0*s50 + 45.0*s10) / 32.0

                xinit = [a, b, c]

            elif self._volatilityFunctionType == FinVolFunctionTypes.CLARK5:

                # Our preference is to fit to the 10D wings first
                if self._useMS10DVol is False:
                    # Fit to 25D
                    c0 = np.log(atmVol)
                    c1 = 2.0 * np.log(s75/s25)
                    c2 = 8.0 * np.log(s25*s75/atmVol/atmVol)
                    xinit = [c0, c1, c2, 0.0, 0.0]
                else:
                    # Fit to 10D
                    c0 = np.log(atmVol)
                    c1 = np.log(s90/s10) / 0.80
                    c2 = np.log(s10*s90/atmVol/atmVol) / 0.32
                    xinit = [c0, c1, c2, 0.0, 0.0]

            else:
                raise FinError("Unknown Model Type")

            xinits.append(xinit)
            ginits.append(ginit)

        deltaMethodValue = self._deltaMethod.value
        volTypeValue = self._volatilityFunctionType.value

        for i in range(0, numVolCurves):

            t = self._texp[i]
            rd = self._rd[i]
            rf = self._rf[i]
            K_ATM = self._K_ATM[i]
            atmVol = self._atmVols[i]

            # If the data has not been provided, pass a dummy value 
            # as I don't want more arguments and Numpy needs floats
            if self._useMS25DVol:
                ms25DVol = self._mktStrangle25DeltaVols[i]
                rr25DVol = self._riskReversal25DeltaVols[i]
            else:
                ms25DVol = -999.0
                rr25DVol = -999.0

            if self._useMS10DVol:
                ms10DVol = self._mktStrangle10DeltaVols[i]
                rr10DVol = self._riskReversal10DeltaVols[i]
            else:
                ms10DVol = -999.0
                rr10DVol = -999.0

            res = _solveToHorizon(s, t, rd, rf,
                                      K_ATM, atmVol,
                                      ms25DVol, rr25DVol,
                                      ms10DVol, rr10DVol,
                                      deltaMethodValue, volTypeValue,
                                      self._alpha,
                                      xinits[i],
                                      ginits[i],
                                      finSolverType,
                                      tol)

            (self._parameters[i,:], self._strikes[i,:], self._gaps[i:],
             self._K_25D_C_MS[i], self._K_25D_P_MS[i],
             self._K_25D_C[i], self._K_25D_P[i],
             self._K_10D_C_MS[i], self._K_10D_P_MS[i],
             self._K_10D_C[i], self._K_10D_P[i]
             ) = res

###############################################################################

    def checkCalibration(self, verbose: bool, tol: float = 1e-6):
        ''' Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals. '''

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valueDate)
            print("SPOT FX RATE:", self._spotFXRate)
            print("ALPHA WEIGHT:", self._alpha)
            print("ATM METHOD:", self._atmMethod)
            print("DELTA METHOD:", self._deltaMethod)
            print("==========================================================")

        K_dummy = 999

        for i in range(0, self._numVolCurves):

            expiryDate = self._expiryDates[i]

            if verbose:
                print("TENOR:", self._tenors[i])
                print("EXPIRY DATE:", expiryDate)
                print("IN ATM VOL: %9.6f %%"%
                      (100.0*self._atmVols[i]))

                if self._useMS25DVol:
                    print("IN MKT STRANGLE 25D VOL: %9.6f %%"%
                          (100.0*self._mktStrangle25DeltaVols[i]))
                    print("IN RSK REVERSAL 25D VOL: %9.6f %%"%
                          (100.0*self._riskReversal25DeltaVols[i]))

                if self._useMS10DVol:
                    print("IN MKT STRANGLE 10D VOL: %9.6f %%"%
                          (100.0*self._mktStrangle10DeltaVols[i]))
                    print("IN RSK REVERSAL 10D VOL: %9.6f %%"%
                          (100.0*self._riskReversal10DeltaVols[i]))

            call = FinFXVanillaOption(expiryDate,
                                      K_dummy,
                                      self._currencyPair,
                                      FinOptionTypes.EUROPEAN_CALL,
                                      1.0,
                                      self._notionalCurrency, )

            put = FinFXVanillaOption(expiryDate,
                                     K_dummy,
                                     self._currencyPair,
                                     FinOptionTypes.EUROPEAN_PUT,
                                     1.0,
                                     self._notionalCurrency)

            ###################################################################
            # AT THE MONEY
            ###################################################################

            if verbose:
                print("==========================================================")
                print("T_(YEARS): ", self._texp[i])
                print("CNT_CPD_RD:%9.6f %%"% (self._rd[i]*100))
                print("CNT_CPD_RF:%9.6f %%"% (self._rf[i]*100))
                print("FWD_RATE:  %9.6f"% (self._F0T[i]))

            sigma_ATM_out = volFunction(self._volatilityFunctionType.value,
                                            self._parameters[i],
                                            self._strikes[i],
                                            self._gaps[i],
                                            self._F0T[i],
                                            self._K_ATM[i],
                                            self._texp[i])

            if verbose:
                print("==========================================================")
                print("VOL FUNCTION", self._volatilityFunctionType)
                print("VOL_PARAMETERS:", self._parameters[i])
                print("==========================================================")
                print("OUT_K_ATM:  %9.6f" % (self._K_ATM[i]))
                print("OUT_ATM_VOL: %9.6f %%"
                      % (100.0*sigma_ATM_out))

            diff = sigma_ATM_out - self._atmVols[i]

            if np.abs(diff) > tol:
                print("FAILED FIT TO ATM VOL IN: %9.6f  OUT: %9.6f  DIFF: %9.6f"%
                      (self._atmVols[i]*100.0, sigma_ATM_out*100.0,
                       diff * 100.0))

            call._strikeFXRate = self._K_ATM[i]
            put._strikeFXRate = self._K_ATM[i]

            model = FinModelBlackScholes(sigma_ATM_out)

            delta_call = call.delta(self._valueDate,
                                    self._spotFXRate,
                                    self._domDiscountCurve,
                                    self._forDiscountCurve,
                                    model)[self._deltaMethodString]

            delta_put = put.delta(self._valueDate,
                                  self._spotFXRate,
                                  self._domDiscountCurve,
                                  self._forDiscountCurve,
                                  model)[self._deltaMethodString]

            if verbose:
                print("CALL_DELTA: % 9.6f  PUT_DELTA: % 9.6f  NET_DELTA: % 9.6f"
                      % (delta_call, delta_put, delta_call + delta_put))

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.25/-0.25
            ###################################################################

            if self._useMS25DVol is True:

                msVol = self._atmVols[i] + self._mktStrangle25DeltaVols[i]

                if verbose:
    
                    print("==========================================================")
                    print("MKT STRANGLE 25D VOL IN: %9.6f %%"
                          % (100.0*self._mktStrangle25DeltaVols[i]))
    
                call._strikeFXRate = self._K_25D_C_MS[i]
                put._strikeFXRate = self._K_25D_P_MS[i]
    
                model = FinModelBlackScholes(msVol)
    
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                if verbose:
                    print("K_25D_C_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_C_MS[i], 100.0*msVol, delta_call))
    
                    print("K_25D_P_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_P_MS[i], 100.0*msVol, delta_put))
    
                call_value = call.value(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)['v']
    
                put_value = put.value(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)['v']
    
                mktStrangleValue = call_value + put_value
    
                if verbose:
                    print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_VALUE: % 9.6f"
                          % (call_value, put_value, mktStrangleValue))
    
                ###################################################################
                # NOW WE ASSIGN A DIFFERENT VOLATILITY TO THE MS STRIKES
                # THE DELTAS WILL NO LONGER EQUAL 0.25, -0.25
                ###################################################################
    
                # CALL
                sigma_K_25D_C_MS = volFunction(self._volatilityFunctionType.value,
                                                   self._parameters[i],
                                                   self._strikes[i],
                                                   self._gaps[i],
                                                   self._F0T[i],
                                                   self._K_25D_C_MS[i],
                                                   self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_25D_C_MS)
                call_value = call.value(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)['v']
    
                # THIS IS NOT GOING TO BE 0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                # PUT
                sigma_K_25D_P_MS = volFunction(self._volatilityFunctionType.value,
                                                   self._parameters[i],
                                                   self._strikes[i],
                                                   self._gaps[i],
                                                   self._F0T[i],
                                                   self._K_25D_P_MS[i],
                                                   self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_25D_P_MS)
                put_value = put.value(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)['v']
    
                # THIS IS NOT GOING TO BE -0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                mktStrangleValueSkew = call_value + put_value
    
                if verbose:
                    print("K_25D_C_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_C_MS[i], 100.0*sigma_K_25D_C_MS, delta_call))
    
                    print("K_25D_P_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_P_MS[i], 100.0*sigma_K_25D_P_MS, delta_put))
    
                    print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_SKEW_VALUE: % 9.6f"
                          % (call_value, put_value, mktStrangleValueSkew))
    
                diff = mktStrangleValue - mktStrangleValueSkew
                if np.abs(diff) > tol:
                    print("FAILED FIT TO 25D MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f"%
                          (mktStrangleValue, mktStrangleValueSkew, diff))
    
                ###################################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.25, -0.25
                ###################################################################
    
                call._strikeFXRate = self._K_25D_C[i]
                put._strikeFXRate = self._K_25D_P[i]
    
                sigma_K_25D_C = volFunction(self._volatilityFunctionType.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_25D_C[i],
                                                self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_25D_C)
    
                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                sigma_K_25D_P = volFunction(self._volatilityFunctionType.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_25D_P[i],
                                                self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_25D_P)
    
                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                if verbose:
                    print("K_25D_C: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                          % (self._K_25D_C[i], 100.0*sigma_K_25D_C, delta_call))
    
                    print("K_25D_P: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                          % (self._K_25D_P[i], 100.0*sigma_K_25D_P, delta_put))
    
                sigma_RR = sigma_K_25D_C - sigma_K_25D_P
    
                if verbose:
                    print("==========================================================")
                    print("RR = VOL_K_25_C - VOL_K_25_P => RR_IN: %9.6f %% RR_OUT: %9.6f %%"
                      % (100.0 * self._riskReversal25DeltaVols[i], 100.0*sigma_RR))
                    print("==========================================================")
    
                diff = sigma_RR - self._riskReversal25DeltaVols[i]
    
                if np.abs(diff) > tol:
                    print("FAILED FIT TO 25D RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f"%
                          (self._riskReversal25DeltaVols[i]*100.0,
                           sigma_RR*100.0,
                           diff*100.0))

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.10/-0.10
            ###################################################################

            if self._useMS10DVol:

                msVol = self._atmVols[i] + self._mktStrangle10DeltaVols[i]
    
                if verbose:
    
                    print("==========================================================")
                    print("MKT STRANGLE 10D VOL IN: %9.6f %%"
                          % (100.0*self._mktStrangle10DeltaVols[i]))
    
                call._strikeFXRate = self._K_10D_C_MS[i]
                put._strikeFXRate = self._K_10D_P_MS[i]
    
                model = FinModelBlackScholes(msVol)
    
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                if verbose:
                    print("K_10D_C_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_C_MS[i], 100.0*msVol, delta_call))
    
                    print("K_10D_P_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_P_MS[i], 100.0*msVol, delta_put))
    
                call_value = call.value(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)['v']
    
                put_value = put.value(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)['v']
    
                mktStrangleValue = call_value + put_value
    
                if verbose:
                    print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_VALUE: % 9.6f"
                          % (call_value, put_value, mktStrangleValue))
    
                ###################################################################
                # NOW WE ASSIGN A DIFFERENT VOLATILITY TO THE MS STRIKES
                # THE DELTAS WILL NO LONGER EQUAL 0.25, -0.25
                ###################################################################
    
                # CALL
                sigma_K_10D_C_MS = volFunction(self._volatilityFunctionType.value,
                                                   self._parameters[i],
                                                   self._strikes[i],
                                                   self._gaps[i],
                                                   self._F0T[i],
                                                   self._K_10D_C_MS[i],
                                                   self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_10D_C_MS)
                call_value = call.value(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)['v']
    
                # THIS IS NOT GOING TO BE 0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                # PUT
                sigma_K_10D_P_MS = volFunction(self._volatilityFunctionType.value,
                                                   self._parameters[i],
                                                   self._strikes[i],
                                                   self._gaps[i],
                                                   self._F0T[i],
                                                   self._K_10D_P_MS[i],
                                                   self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_10D_P_MS)
                put_value = put.value(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)['v']
    
                # THIS IS NOT GOING TO BE -0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                mktStrangleValueSkew = call_value + put_value
    
                if verbose:
                    print("K_10D_C_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_C_MS[i], 100.0*sigma_K_10D_C_MS, delta_call))
    
                    print("K_10D_P_MS: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_P_MS[i], 100.0*sigma_K_10D_P_MS, delta_put))
    
                    print("CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_SKEW_VALUE: % 9.6f"
                          % (call_value, put_value, mktStrangleValueSkew))
    
                diff = mktStrangleValue - mktStrangleValueSkew
                if np.abs(diff) > tol:
                    print("FAILED FIT TO 10D MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f"%
                          (mktStrangleValue, mktStrangleValueSkew, diff))
    
                ###################################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.10, -0.10
                ###################################################################
    
                call._strikeFXRate = self._K_10D_C[i]
                put._strikeFXRate = self._K_10D_P[i]
    
                sigma_K_10D_C = volFunction(self._volatilityFunctionType.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_10D_C[i],
                                                self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_10D_C)
    
                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(self._valueDate,
                                        self._spotFXRate,
                                        self._domDiscountCurve,
                                        self._forDiscountCurve,
                                        model)[self._deltaMethodString]
    
                sigma_K_10D_P = volFunction(self._volatilityFunctionType.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_10D_P[i],
                                                self._texp[i])
    
                model = FinModelBlackScholes(sigma_K_10D_P)
    
                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(self._valueDate,
                                      self._spotFXRate,
                                      self._domDiscountCurve,
                                      self._forDiscountCurve,
                                      model)[self._deltaMethodString]
    
                if verbose:
                    print("K_10D_C: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                          % (self._K_10D_C[i], 100.0*sigma_K_10D_C, delta_call))
    
                    print("K_10D_P: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                          % (self._K_10D_P[i], 100.0*sigma_K_10D_P, delta_put))
    
                sigma_RR = sigma_K_10D_C - sigma_K_10D_P
    
                if verbose:
                    print("==========================================================")
                    print("RR = VOL_K_10D_C - VOL_K_10D_P => RR_IN: %9.6f %% RR_OUT: %9.6f %%"
                      % (100.0 * self._riskReversal10DeltaVols[i], 100.0*sigma_RR))
                    print("==========================================================")
    
                diff = sigma_RR - self._riskReversal10DeltaVols[i]
    
                if np.abs(diff) > tol:
                    print("FAILED FIT TO 10D RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f"%
                          (self._riskReversal10DeltaVols[i]*100.0,
                           sigma_RR*100.0,
                           diff*100.0))
    
###############################################################################

    def impliedDbns(self, lowFX, highFX, numIntervals):
        ''' Calculate the pdf for each tenor horizon. Returns a list of
        FinDistribution objects, one for each tenor horizon. '''

        dbns = []

        for iTenor in range(0, len(self._tenors)):

            f = self._F0T[iTenor]
            t = self._texp[iTenor]

            dFX = (highFX - lowFX)/ numIntervals

            domDF = self._domDiscountCurve._df(t)
            forDF = self._forDiscountCurve._df(t)

            rd = -np.log(domDF) / t
            rf = -np.log(forDF) / t

            Ks = []
            vols = []

            for iK in range(0, numIntervals):

                k = lowFX + iK*dFX

                vol = volFunction(self._volatilityFunctionType.value, 
                                      self._parameters[iTenor], 
                                      self._strikes[iTenor], 
                                      self._gaps[iTenor], 
                                      f, k, t)

                Ks.append(k)
                vols.append(vol)

            Ks = np.array(Ks)
            vols = np.array(vols)

            density = optionImpliedDbn(self._spotFXRate, t, rd, rf, Ks, vols)

            dbn = FinDistribution(Ks, density)
            dbns.append(dbn)

        return dbns

###############################################################################

    def plotVolCurves(self):
        ''' Generates a plot of each of the vol curves implied by the market 
        and fitted. '''
        
        plt.figure()

        volTypeVal = self._volatilityFunctionType.value

        for tenorIndex in range(0, self._numVolCurves):

            atmVol = self._atmVols[tenorIndex]*100
            msVol25 = self._mktStrangle25DeltaVols[tenorIndex]*100
            rrVol25 = self._riskReversal25DeltaVols[tenorIndex]*100
            msVol10 = self._mktStrangle10DeltaVols[tenorIndex]*100
            rrVol10 = self._riskReversal10DeltaVols[tenorIndex]*100
            strikes = self._strikes[tenorIndex]
            
            gaps = self._gaps[tenorIndex]

            lowK = self._K_10D_P[tenorIndex] * 0.90
            highK = self._K_10D_C_MS[tenorIndex] * 1.10

            ks = []
            vols = []
            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            params = self._parameters[tenorIndex]
            t = self._texp[tenorIndex]
            f = self._F0T[tenorIndex]
            
            for i in range(0, numIntervals):

                sigma = volFunction(volTypeVal, params, strikes, gaps,
                                    f, K, t) * 100.0
                ks.append(K)
                vols.append(sigma)
                K = K + dK

            labelStr = self._tenors[tenorIndex]
            labelStr += " ATM: " + str(atmVol)[0:6]
            labelStr += " MS25: " + str(msVol25)[0:6]
            labelStr += " RR25: " + str(rrVol25)[0:6]
            labelStr += " MS10: " + str(msVol10)[0:6]
            labelStr += " RR10: " + str(rrVol10)[0:6]

            plt.plot(ks, vols, label=labelStr)
            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = "JNT FIT:" + self._currencyPair + " " +\
                    str(self._volatilityFunctionType)

            keyStrikes = []
            keyStrikes.append(self._K_ATM[tenorIndex])

            keyVols = []
            for K in keyStrikes:

                sigma = volFunction(volTypeVal, params, 
                                        strikes, gaps,
                                        f, K, t) * 100.0

                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'ko', markersize=4)

            keyStrikes = []
            keyStrikes.append(self._K_25D_P[tenorIndex])
            keyStrikes.append(self._K_25D_P_MS[tenorIndex])
            keyStrikes.append(self._K_25D_C[tenorIndex])
            keyStrikes.append(self._K_25D_C_MS[tenorIndex])

            keyVols = []
            for K in keyStrikes:

                sigma = volFunction(volTypeVal, params, 
                                        strikes, gaps,
                                        f, K, t) * 100.0

                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'bo', markersize=4)

            keyStrikes = []
            keyStrikes.append(self._K_10D_P[tenorIndex])
            keyStrikes.append(self._K_10D_P_MS[tenorIndex])
            keyStrikes.append(self._K_10D_C[tenorIndex])
            keyStrikes.append(self._K_10D_C_MS[tenorIndex])

            keyVols = []
            for K in keyStrikes:
                sigma = volFunction(volTypeVal, params, 
                                        strikes, gaps, 
                                        f, K, t) * 100.0
                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'ro', markersize=4)

        plt.title(title)
        plt.legend(loc="lower left", bbox_to_anchor=(1,0))

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VALUE DATE", self._valueDate)
        s += labelToString("FX RATE", self._spotFXRate)
        s += labelToString("CCY PAIR", self._currencyPair)
        s += labelToString("NOTIONAL CCY", self._notionalCurrency)
        s += labelToString("NUM TENORS", self._numVolCurves)
        s += labelToString("ATM METHOD", self._atmMethod)
        s += labelToString("DELTA METHOD", self._deltaMethod)
        s += labelToString("ALPHA WEIGHT", self._alpha)
        s += labelToString("VOL FUNCTION", self._volatilityFunctionType)

        for i in range(0, self._numVolCurves):

            s += "\n"

            s += labelToString("TENOR", self._tenors[i])
            s += labelToString("EXPIRY DATE", self._expiryDates[i])
            s += labelToString("TIME (YRS)", self._texp[i])
            s += labelToString("FWD FX", self._F0T[i])

            s += labelToString("ATM VOLS", self._atmVols[i]*100.0)
            s += labelToString("MS VOLS", self._mktStrangle25DeltaVols[i]*100.)
            s += labelToString("RR VOLS", self._riskReversal25DeltaVols[i]*100.)

            s += labelToString("ATM Strike", self._K_ATM[i])
            s += labelToString("ATM Delta", self._deltaATM[i])

            s += labelToString("K_ATM", self._K_ATM[i])

            s += labelToString("MS 25D Call Strike", self._K_25D_C_MS[i])
            s += labelToString("MS 25D Put Strike", self._K_25D_P_MS[i])
            s += labelToString("SKEW 25D CALL STRIKE", self._K_25D_C[i])
            s += labelToString("SKEW 25D PUT STRIKE", self._K_25D_P[i])
            s += labelToString("PARAMS", self._parameters[i])

            s += labelToString("MS 10D Call Strike", self._K_10D_C_MS[i])
            s += labelToString("MS 10D Put Strike", self._K_10D_P_MS[i])
            s += labelToString("SKEW 10D CALL STRIKE", self._K_10D_C[i])
            s += labelToString("SKEW 10D PUT STRIKE", self._K_10D_P[i])

        return s

###############################################################################

    def _print(self):
        ''' Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. '''
        print(self)

###############################################################################
