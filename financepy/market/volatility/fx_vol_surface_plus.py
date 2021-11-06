##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import gDaysInYear
from ...utils.global_types import OptionTypes
from ...products.fx.fx_vanilla_option import FXVanillaOption
from ...models.option_implied_dbn import option_implied_dbn
from ...products.fx.fx_mkt_conventions import FinFXATMMethod
from ...products.fx.fx_mkt_conventions import FinFXDeltaMethod
from ...utils.helpers import check_argument_types, label_to_string
from ...market.curves.discount_curve import DiscountCurve

from ...models.black_scholes import BlackScholes

from ...models.volatility_fns import vol_function_clark
from ...models.volatility_fns import vol_function_bloomberg
from ...models.volatility_fns import VolFunctionTypes
from ...models.sabr import vol_function_sabr
from ...models.sabr import vol_function_sabr_beta_one
from ...models.sabr import vol_function_sabr_beta_half

from ...utils.math import norminvcdf

from ...models.black_scholes_analytic import bs_value
from ...products.fx.fx_vanilla_option import fast_delta
from ...utils.distribution import FinDistribution

from ...utils.solver_1d import newton_secant
from ...utils.solver_nm import nelder_mead
from ...utils.global_types import FinSolverTypes

###############################################################################
# ISSUES
# sabr does not fit inverted skew discount like eurjpy
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
    """ This is the objective function used in the determination of the FX
    option implied strike which is computed in the class below. """

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    volatility = args[4]
    delta_method_value = args[5]
    option_type_value = args[6]
    delta_target = args[7]

    delta_out = fast_delta(s, t, K, rd, rf,
                           volatility,
                           delta_method_value,
                           option_type_value)

    obj_fn = delta_target - delta_out
    return obj_fn

###############################################################################


@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def _interpolate_gap(k, strikes, gaps):

    if k <= strikes[0]:
        return 0.0

    if k >= strikes[-1]:
        return 0.0

    index = 0
    for i in range(1, len(strikes)):
        if k > strikes[i-1] and k <= strikes[i]:
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


@njit(fastmath=True)  # , cache=True)
def _obj(params, *args):
    """ Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the vol_type_value
    """

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    K_ATM = args[4]
    atm_vol = args[5]

    K_25D_C_MS = args[6]
    K_25D_P_MS = args[7]
    V_25D_MS_target = args[8]
    target25DRRVol = args[9]

    K_10D_C_MS = args[10]
    K_10D_P_MS = args[11]
    V_10D_MS_target = args[12]
    target10DRRVol = args[13]

    delta_method_value = args[14]
    vol_type_value = args[15]
    alpha = args[16]

    strikesNULL = np.zeros(1)
    gapsNULL = np.zeros(1)

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atm_curve_vol = vol_function(vol_type_value, params, strikesNULL, gapsNULL,
                                 f, K_ATM, t)

    termATM = (atm_vol - atm_curve_vol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25D strikes
    ###########################################################################

    if target25DRRVol > -999.0:

        sigma_K_25D_C_MS = vol_function(vol_type_value, params,
                                        strikesNULL, gapsNULL,
                                        f, K_25D_C_MS, t)

        V_25D_C_MS = bs_value(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                              OptionTypes.EUROPEAN_CALL.value)

        sigma_K_25D_P_MS = vol_function(vol_type_value, params,
                                        strikesNULL, gapsNULL,
                                        f, K_25D_P_MS, t)

        V_25D_P_MS = bs_value(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                              OptionTypes.EUROPEAN_PUT.value)

        V_25D_MS = V_25D_C_MS + V_25D_P_MS
        term25D_1 = (V_25D_MS - V_25D_MS_target)**2

    else:

        term25D_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target25DRRVol > -999.0:

        K_25D_C = _solver_for_smile_strike(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_CALL.value,
                                           vol_type_value, +0.2500,
                                           delta_method_value, K_25D_C_MS,
                                           params, strikesNULL, gapsNULL)

        sigma_K_25D_C = vol_function(vol_type_value, params,
                                     strikesNULL, gapsNULL,
                                     f, K_25D_C, t)

        K_25D_P = _solver_for_smile_strike(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_PUT.value,
                                           vol_type_value, -0.2500,
                                           delta_method_value, K_25D_P_MS,
                                           params, strikesNULL, gapsNULL)

        sigma_K_25D_P = vol_function(vol_type_value, params,
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

        sigma_K_10D_C_MS = vol_function(vol_type_value, params,
                                        strikesNULL, gapsNULL,
                                        f, K_10D_C_MS, t)

        V_10D_C_MS = bs_value(s, t, K_10D_C_MS, rd, rf, sigma_K_10D_C_MS,
                              OptionTypes.EUROPEAN_CALL.value)

        sigma_K_10D_P_MS = vol_function(vol_type_value, params,
                                        strikesNULL, gapsNULL,
                                        f, K_10D_P_MS, t)

        V_10D_P_MS = bs_value(s, t, K_10D_P_MS, rd, rf, sigma_K_10D_P_MS,
                              OptionTypes.EUROPEAN_PUT.value)

        V_10D_MS = V_10D_C_MS + V_10D_P_MS
        term10D_1 = (V_10D_MS - V_10D_MS_target)**2

    else:

        term10D_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target10DRRVol > -999.0:

        K_10D_C = _solver_for_smile_strike(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_CALL.value,
                                           vol_type_value, +0.1000,
                                           delta_method_value, K_10D_C_MS,
                                           params, strikesNULL, gapsNULL)

        sigma_K_10D_C = vol_function(vol_type_value, params,
                                     strikesNULL, gapsNULL,
                                     f, K_10D_C, t)

        K_10D_P = _solver_for_smile_strike(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_PUT.value,
                                           vol_type_value, -0.1000,
                                           delta_method_value, K_10D_P_MS,
                                           params, strikesNULL, gapsNULL)

        sigma_K_10D_P = vol_function(vol_type_value, params,
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


@njit(fastmath=True)  # , cache=True)
def _obj_gap(gaps, *args):
    """ Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the vol_type_value
    """

    s = args[0]
    t = args[1]
    rd = args[2]
    rf = args[3]
    K_ATM = args[4]
    atm_vol = args[5]

    K_25D_C_MS = args[6]
    K_25D_P_MS = args[7]
    V_25D_MS_target = args[8]
    target25DRRVol = args[9]

    K_10D_C_MS = args[10]
    K_10D_P_MS = args[11]
    V_10D_MS_target = args[12]
    target10DRRVol = args[13]

    delta_method_value = args[14]
    vol_type_value = args[15]
    params = args[16]

    strikes = [K_10D_P_MS, K_25D_P_MS, K_ATM, K_25D_C_MS, K_10D_C_MS]
    strikes = np.array(strikes)

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atm_curve_vol = vol_function(vol_type_value, params, strikes, gaps,
                                 f, K_ATM, t)

    print("atm_curve_vol", atm_curve_vol)

    termATM = (atm_vol - atm_curve_vol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25D strikes
    ###########################################################################

    sigma_K_25D_C_MS = vol_function(vol_type_value, params, strikes, gaps,
                                    f, K_25D_C_MS, t)

    print("sigma_K_25D_C_MS", sigma_K_25D_C_MS)

    V_25D_C_MS = bs_value(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                          OptionTypes.EUROPEAN_CALL.value)

    sigma_K_25D_P_MS = vol_function(vol_type_value, params, strikes, gaps,
                                    f, K_25D_P_MS, t)

    print("sigma_K_25D_P_MS", sigma_K_25D_P_MS)

    V_25D_P_MS = bs_value(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                          OptionTypes.EUROPEAN_PUT.value)

    V_25D_MS = V_25D_C_MS + V_25D_P_MS
    term25D_1 = (V_25D_MS - V_25D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    K_25D_C = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_CALL.value,
                                       vol_type_value, +0.2500,
                                       delta_method_value, K_25D_C_MS,
                                       params, strikes, gaps)

    sigma_K_25D_C = vol_function(vol_type_value, params, strikes, gaps,
                                 f, K_25D_C, t)

    print("sigma_K_25D_C", sigma_K_25D_C)

    K_25D_P = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_PUT.value,
                                       vol_type_value, -0.2500,
                                       delta_method_value, K_25D_P_MS,
                                       params, strikes, gaps)

    sigma_K_25D_P = vol_function(vol_type_value, params, strikes, gaps,
                                 f, K_25D_P, t)

    print("sigma_K_25D_P", sigma_K_25D_P)

    sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
    term25D_2 = (sigma_25D_RR - target25DRRVol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 10D strikes
    ###########################################################################

    sigma_K_10D_C_MS = vol_function(vol_type_value, params, strikes, gaps,
                                    f, K_10D_C_MS, t)

    print("sigma_K_10D_C_MS", sigma_K_10D_C_MS)

    V_10D_C_MS = bs_value(s, t, K_10D_C_MS, rd, rf, sigma_K_10D_C_MS,
                          OptionTypes.EUROPEAN_CALL.value)

    sigma_K_10D_P_MS = vol_function(vol_type_value, params, strikes, gaps,
                                    f, K_10D_P_MS, t)

    print("sigma_K_10D_P_MS", sigma_K_10D_P_MS)

    V_10D_P_MS = bs_value(s, t, K_10D_P_MS, rd, rf, sigma_K_10D_P_MS,
                          OptionTypes.EUROPEAN_PUT.value)

    V_10D_MS = V_10D_C_MS + V_10D_P_MS
    term10D_1 = (V_10D_MS - V_10D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    K_10D_C = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_CALL.value,
                                       vol_type_value, +0.1000,
                                       delta_method_value, K_10D_C_MS,
                                       params, strikes, gaps)

    sigma_K_10D_C = vol_function(vol_type_value, params, strikes, gaps,
                                 f, K_10D_C, t)

    print("SIGMA_K_10D_C", sigma_K_10D_C)

    print("INIT K_10D_P_MS", K_10D_P_MS)

    K_10D_P = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_PUT.value,
                                       vol_type_value, -0.1000,
                                       delta_method_value, K_10D_P_MS,
                                       params, strikes, gaps)

    print("K_10D_P", K_10D_P)
    sigma_K_10D_P = vol_function(vol_type_value, params, strikes, gaps,
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


def _solve_to_horizon(s, t, rd, rf,
                      K_ATM, atm_vol,
                      ms25DVol, rr25DVol,
                      ms10DVol, rr10DVol,
                      delta_method_value, vol_type_value,
                      alpha,
                      x_inits,
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

        vol_25D_MS = atm_vol + ms25DVol

        K_25D_C_MS = solve_for_strike(s, t, rd, rf,
                                      OptionTypes.EUROPEAN_CALL.value,
                                      +0.2500,
                                      delta_method_value,
                                      vol_25D_MS)

        K_25D_P_MS = solve_for_strike(s, t, rd, rf,
                                      OptionTypes.EUROPEAN_PUT.value,
                                      -0.2500,
                                      delta_method_value,
                                      vol_25D_MS)

        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        V_25D_C_MS = bs_value(s, t, K_25D_C_MS, rd, rf, vol_25D_MS,
                              OptionTypes.EUROPEAN_CALL.value)

        V_25D_P_MS = bs_value(s, t, K_25D_P_MS, rd, rf, vol_25D_MS,
                              OptionTypes.EUROPEAN_PUT.value)

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

        vol_10D_MS = atm_vol + ms10DVol

        K_10D_C_MS = solve_for_strike(s, t, rd, rf,
                                      OptionTypes.EUROPEAN_CALL.value,
                                      +0.1000,
                                      delta_method_value,
                                      vol_10D_MS)

        K_10D_P_MS = solve_for_strike(s, t, rd, rf,
                                      OptionTypes.EUROPEAN_PUT.value,
                                      -0.1000,
                                      delta_method_value,
                                      vol_10D_MS)

        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        V_10D_C_MS = bs_value(s, t, K_10D_C_MS, rd, rf, vol_10D_MS,
                              OptionTypes.EUROPEAN_CALL.value)

        V_10D_P_MS = bs_value(s, t, K_10D_P_MS, rd, rf, vol_10D_MS,
                              OptionTypes.EUROPEAN_PUT.value)

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
            K_ATM, atm_vol,
            K_25D_C_MS, K_25D_P_MS, V_25D_MS, rr25DVol,
            K_10D_C_MS, K_10D_P_MS, V_10D_MS, rr10DVol,
            delta_method_value, vol_type_value, alpha)

    # Nelmer-Mead (both SciPy & Numba) is quicker, but occasionally fails
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if finSolverType == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(_obj, np.array(x_inits),
                               bounds=np.array(
                                   [[], []]).T, args=args, tol_f=tol,
                               tol_x=tol, max_iter=1000)
        elif finSolverType == FinSolverTypes.NELDER_MEAD:
            opt = minimize(_obj, x_inits, args, method="Nelder-Mead", tol=tol)
            xopt = opt.x
        elif finSolverType == FinSolverTypes.CONJUGATE_GRADIENT:
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)
            xopt = opt.x
    except:
        # If convergence fails try again with CG if necessary
        if finSolverType != FinSolverTypes.CONJUGATE_GRADIENT:
            print('Failed to converge, will try CG')
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)

            xopt = opt.x

    params = np.array(xopt)

    strikes = [K_10D_P_MS, K_25D_P_MS, K_ATM, K_10D_C_MS, K_25D_C_MS]
    strikes = np.array(strikes)
    gaps = np.zeros(5)

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    if 1 == 0:

        tol = 1e-12

        args = (s, t, rd, rf,
                K_ATM, atm_vol,
                K_25D_C_MS, K_25D_P_MS, V_25D_MS, rr25DVol,
                K_10D_C_MS, K_10D_P_MS, V_10D_MS, rr10DVol,
                delta_method_value, vol_type_value, params)

        opt = minimize(_obj_gap, ginits, args, method="Nelder-Mead", tol=tol)
        xopt = opt.x
        gaps = np.array(xopt)

        print("SOLVED")

# Removed this as it causes discontinuity
#    f = s * np.exp((rd-rf)*t)
#    interpATMVol = vol_function(vol_type_value, params,
#                                   strikes, gaps, f, K_ATM, t)

#    diff = atm_vol - interpATMVol
#    gaps[2] = diff

    ###########################################################################

    if use25D is False:
        K_25D_C_MS = K_ATM
        K_25D_P_MS = K_ATM

    K_25D_C = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_CALL.value,
                                       vol_type_value, +0.2500,
                                       delta_method_value, K_25D_C_MS,
                                       params, strikes, gaps)

    K_25D_P = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_PUT.value,
                                       vol_type_value, -0.2500,
                                       delta_method_value, K_25D_P_MS,
                                       params, strikes, gaps)

    if use10D is False:
        K_10D_C_MS = K_ATM
        K_10D_P_MS = K_ATM

    K_10D_C = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_CALL.value,
                                       vol_type_value, +0.1000,
                                       delta_method_value, K_10D_C_MS,
                                       params, strikes, gaps)

    K_10D_P = _solver_for_smile_strike(s, t, rd, rf,
                                       OptionTypes.EUROPEAN_PUT.value,
                                       vol_type_value, -0.1000,
                                       delta_method_value, K_10D_P_MS,
                                       params, strikes, gaps)

    return (params, strikes, gaps,
            K_25D_C_MS, K_25D_P_MS, K_25D_C, K_25D_P,
            K_10D_C_MS, K_10D_P_MS, K_10D_C, K_10D_P)

###############################################################################


@njit(float64(int64, float64[:], float64[:], float64[:],
              float64, float64, float64), cache=True, fastmath=True)
def vol_function(vol_function_type_value, params, strikes, gaps, f, k, t):
    """ Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. """

#    print("vol_function", vol_function_type_value)

    if len(strikes) == 1:
        gapK = 0.0
    else:
        gapK = _interpolate_gap(k, strikes, gaps)

    if vol_function_type_value == VolFunctionTypes.CLARK.value:
        vol = vol_function_clark(params, f, k, t) + gapK
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR.value:
        vol = vol_function_sabr(params, f, k, t) + gapK
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_HALF.value:
        vol = vol_function_sabr_beta_half(params, f, k, t) + gapK
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_ONE.value:
        vol = vol_function_sabr_beta_one(params, f, k, t) + gapK
        return vol
    elif vol_function_type_value == VolFunctionTypes.BBG.value:
        vol = vol_function_bloomberg(params, f, k, t) + gapK
        return vol
    elif vol_function_type_value == VolFunctionTypes.CLARK5.value:
        vol = vol_function_clark(params, f, k, t) + gapK
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


@njit(cache=True, fastmath=True)
def _delta_fit(k, *args):
    """ This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol
    and high K. It speeds up the code as it allows initial values close to
    the solution to be used. """

    vol_type_value = args[0]
    s = args[1]
    t = args[2]
    rd = args[3]
    rf = args[4]
    option_type_value = args[5]
    deltaTypeValue = args[6]
    inverseDeltaTarget = args[7]
    params = args[8]
    strikes = args[9]
    gaps = args[10]

    f = s * np.exp((rd-rf)*t)
    v = vol_function(vol_type_value, params, strikes, gaps, f, k, t)
    delta_out = fast_delta(
        s, t, k, rd, rf, v, deltaTypeValue, option_type_value)
    inverseDeltaOut = norminvcdf(np.abs(delta_out))
    invObjFn = inverseDeltaTarget - inverseDeltaOut

#    print(k, f, v, delta_out, invObjFn)

    return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, int64, float64,
              int64, float64, float64[:], float64[:], float64[:]),
      fastmath=True)
def _solver_for_smile_strike(s, t, rd, rf,
                             option_type_value,
                             volatilityTypeValue,
                             delta_target,
                             delta_method_value,
                             initialGuess,
                             parameters,
                             strikes,
                             gaps):
    """ Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike. """

    inverseDeltaTarget = norminvcdf(np.abs(delta_target))

    argtuple = (volatilityTypeValue, s, t, rd, rf,
                option_type_value, delta_method_value,
                inverseDeltaTarget,
                parameters, strikes, gaps)

    K = newton_secant(_delta_fit, x0=initialGuess, args=argtuple,
                      tol=1e-8, maxiter=50)

    return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


@njit(float64(float64, float64, float64, float64, int64, float64,
              int64, float64), fastmath=True)
def solve_for_strike(spot_fx_rate,
                     tdel, rd, rf,
                     option_type_value,
                     delta_target,
                     delta_method_value,
                     volatility):
    """ This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorith to determine the strike that matches an
    input volatility. """

    # =========================================================================
    # IMPORTANT NOTE:
    # =========================================================================
    # For some delta quotation conventions I can solve for K explicitly.
    # Note that as I am using the function norminvdelta to calculate the
    # inverse value of delta, this may not, on a round trip using N(x), give
    # back the value x as it is calculated to a different number of decimal
    # places. It should however agree to 6-7 decimal places. Which is OK.
    # =========================================================================

    if delta_method_value == FinFXDeltaMethod.SPOT_DELTA.value:

        domDF = np.exp(-rd*tdel)
        forDF = np.exp(-rf*tdel)

        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spot_fx_rate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = delta_target*phi/forDF  # CHECK THIS !!!
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
        return K

    elif delta_method_value == FinFXDeltaMethod.FORWARD_DELTA.value:

        domDF = np.exp(-rd*tdel)
        forDF = np.exp(-rf*tdel)

        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        F0T = spot_fx_rate * forDF / domDF
        vsqrtt = volatility * np.sqrt(tdel)
        arg = delta_target*phi
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt/2.0))
        return K

    elif delta_method_value == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

        argtuple = (spot_fx_rate, tdel, rd, rf, volatility,
                    delta_method_value, option_type_value, delta_target)

        K = newton_secant(_g, x0=spot_fx_rate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    elif delta_method_value == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

        argtuple = (spot_fx_rate, tdel, rd, rf, volatility,
                    delta_method_value, option_type_value, delta_target)

        K = newton_secant(_g, x0=spot_fx_rate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    else:

        raise FinError("Unknown FinFXDeltaMethod")

###############################################################################


class FXVolSurfacePlus():
    """ Class to perform a calibration of a chosen parametrised surface to the
    prices of FX options at different strikes and expiry tenors. The
    calibration inputs are the ATM and 25 and 10 Delta volatilities in terms of
    the market strangle amd risk reversals. There is a choice of volatility
    function from cubic in delta to full SABR. Check out VolFunctionTypes.
    Parameter alpha [0,1] is used to interpolate between fitting only 25D when
    alpha=0 to fitting only 10D when alpha=1.0. Alpha=0.5 assigns equal weights
    A vol function with more parameters will give a better fit. Of course. But 
    it might also overfit. Visualising the volatility curve is useful. Also, 
    there is no guarantee that the implied pdf will be positive."""

    def __init__(self,
                 valuation_date: Date,
                 spot_fx_rate: float,
                 currency_pair: str,
                 notional_currency: str,
                 dom_discount_curve: DiscountCurve,
                 for_discount_curve: DiscountCurve,
                 tenors: (list),
                 atm_vols: (list, np.ndarray),
                 mktStrangle25DeltaVols: (list, np.ndarray),
                 riskReversal25DeltaVols: (list, np.ndarray),
                 mktStrangle10DeltaVols: (list, np.ndarray),
                 riskReversal10DeltaVols: (list, np.ndarray),
                 alpha: float,
                 atmMethod: FinFXATMMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL,
                 deltaMethod: FinFXDeltaMethod = FinFXDeltaMethod.SPOT_DELTA,
                 volatility_function_type: VolFunctionTypes = VolFunctionTypes.CLARK,
                 finSolverType: FinSolverTypes = FinSolverTypes.NELDER_MEAD,
                 tol: float = 1e-8):
        """ Create the FinFXVolSurfacePlus object by passing in market vol data
        for ATM, 25 Delta and 10 Delta strikes. The alpha weight shifts the
        fitting between 25D and 10D. Alpha = 0.0 is 100% 25D while alpha = 1.0
        is 100% 10D. An alpha of 0.50 is equally weighted. """

        # I want to allow Nones for some of the market inputs
        if mktStrangle10DeltaVols is None:
            mktStrangle10DeltaVols = []

        if riskReversal10DeltaVols is None:
            riskReversal10DeltaVols = []

        if mktStrangle25DeltaVols is None:
            mktStrangle25DeltaVols = []

        if riskReversal25DeltaVols is None:
            riskReversal25DeltaVols = []

        check_argument_types(self.__init__, locals())

        self._valuation_date = valuation_date
        self._spot_fx_rate = spot_fx_rate
        self._currency_pair = currency_pair

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self._forName = self._currency_pair[0:3]
        self._domName = self._currency_pair[3:6]

        self._notional_currency = notional_currency
        self._dom_discount_curve = dom_discount_curve
        self._for_discount_curve = for_discount_curve
        self._num_vol_curves = len(tenors)
        self._tenors = tenors

        if len(atm_vols) != self._num_vol_curves:
            raise FinError("Number ATM vols must equal number of tenors")

        self._atm_vols = np.array(atm_vols)/100.0

        self._useMS25DVol = True
        self._useRR25DVol = True
        self._useMS10DVol = True
        self._useRR10DVol = True

        # Some of these can be missing which is signified by length zero
        n = len(mktStrangle25DeltaVols)

        if n != self._num_vol_curves and n != 0:
            raise FinError("Number MS25D vols must equal number of tenors")

        if n == 0:
            self._useMS25DVol = False

        n = len(riskReversal25DeltaVols)

        if n != self._num_vol_curves and n != 0:
            raise FinError("Number RR25D vols must equal number of tenors")

        if n == 0:
            self._useRR25DVol = False

        n = len(mktStrangle10DeltaVols)

        if n != self._num_vol_curves and n != 0:
            raise FinError("Number MS10D vols must equal number of tenors")

        if n == 0:
            self._useMS10DVol = False

        n = len(riskReversal10DeltaVols)

        if n != self._num_vol_curves and n != 0:
            raise FinError("Number RR10D vols must equal number of tenors")

        if n == 0:
            self._useRR10DVol = False

        if self._useMS10DVol != self._useRR10DVol:
            raise FinError("You must provide both 10D RR + 10D MS or neither")

        if self._useMS25DVol != self._useRR25DVol:
            raise FinError("You must provide both 25D RR + 25D MS or neither")

        if self._useMS10DVol is False and self._useMS25DVol is False:
            raise FinError(
                "No MS and RR. You must provide 10D or 25D MS + RR.")

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

        self._volatility_function_type = volatility_function_type
        self._tenorIndex = 0

        self._expiry_dates = []
        for i in range(0, self._num_vol_curves):
            expiry_date = valuation_date.add_tenor(tenors[i])
            self._expiry_dates.append(expiry_date)

        self._build_vol_surface(finSolverType=finSolverType, tol=tol)

###############################################################################

    def volatility_from_strike_date(self, K, expiry_date):
        """ Interpolates the Black-Scholes volatility from the volatility
        surface given call option strike and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility."""

        texp = (expiry_date - self._valuation_date) / gDaysInYear

        vol_type_value = self._volatility_function_type.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._num_vol_curves

        if num_curves == 1:

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

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]

        t0 = self._texp[index0]
        t1 = self._texp[index1]

        vol0 = vol_function(vol_type_value, self._parameters[index0],
                            self._strikes[index0], self._gaps[index0],
                            fwd0, K, t0)

        if index1 != index0:

            vol1 = vol_function(vol_type_value, self._parameters[index1],
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

    def delta_to_strike(self, callDelta, expiry_date, deltaMethod):
        """ Interpolates the strike at a delta and expiry date. Linear
        interpolation is used in strike."""

        texp = (expiry_date - self._valuation_date) / gDaysInYear

        vol_type_value = self._volatility_function_type.value

        s = self._spot_fx_rate

        if deltaMethod is None:
            delta_method_value = self._deltaMethod.value
        else:
            delta_method_value = deltaMethod.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._num_vol_curves

        # If there is only one time horizon then assume flat vol to this time
        if num_curves == 1:

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

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        #######################################################################

        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._K_ATM[index0]

        K0 = _solver_for_smile_strike(s, texp, self._rd[index0], self._rf[index0],
                                      OptionTypes.EUROPEAN_CALL.value,
                                      vol_type_value, callDelta,
                                      delta_method_value,
                                      initialGuess,
                                      self._parameters[index0],
                                      self._strikes[index0],
                                      self._gaps[index0])

        if index1 != index0:

            K1 = _solver_for_smile_strike(s, texp,
                                          self._rd[index1],
                                          self._rf[index1],
                                          OptionTypes.EUROPEAN_CALL.value,
                                          vol_type_value, callDelta,
                                          delta_method_value,
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

    def volatility_from_delta_date(self, callDelta, expiry_date,
                                   deltaMethod=None):
        """ Interpolates the Black-Scholes volatility from the volatility
        surface given a call option delta and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are 
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be 
        overriden by a provided delta convention. The resulting volatilities 
        are then determined for each bracketing expiry time and linear 
        interpolation is done in variance space and then converted back to a 
        lognormal volatility."""

        texp = (expiry_date - self._valuation_date) / gDaysInYear

        vol_type_value = self._volatility_function_type.value

        s = self._spot_fx_rate

        if deltaMethod is None:
            delta_method_value = self._deltaMethod.value
        else:
            delta_method_value = deltaMethod.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._num_vol_curves

        # If there is only one time horizon then assume flat vol to this time
        if num_curves == 1:

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

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if texp <= self._texp[i] and texp > self._texp[i-1]:
                    index0 = i-1
                    index1 = i
                    break

        fwd0 = self._F0T[index0]
        fwd1 = self._F0T[index1]

        t0 = self._texp[index0]
        t1 = self._texp[index1]

        initialGuess = self._K_ATM[index0]

        K0 = _solver_for_smile_strike(s, texp, self._rd[index0], self._rf[index0],
                                      OptionTypes.EUROPEAN_CALL.value,
                                      vol_type_value, callDelta,
                                      delta_method_value,
                                      initialGuess,
                                      self._parameters[index0],
                                      self._strikes[index0],
                                      self._gaps[index0])

        vol0 = vol_function(vol_type_value, self._parameters[index0],
                            self._strikes[index0], self._gaps[index0],
                            fwd0, K0, t0)

        if index1 != index0:

            K1 = _solver_for_smile_strike(s, texp,
                                          self._rd[index1],
                                          self._rf[index1],
                                          OptionTypes.EUROPEAN_CALL.value,
                                          vol_type_value, callDelta,
                                          delta_method_value,
                                          initialGuess,
                                          self._parameters[index1],
                                          self._strikes[index1],
                                          self._gaps[index1])

            vol1 = vol_function(vol_type_value, self._parameters[index1],
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
                raise FinError(
                    "Failed interpolation due to negative variance.")

            volt = np.sqrt(vart/texp)

        else:

            volt = vol0
            kt = K0

        return volt, kt

###############################################################################

    def _build_vol_surface(self, finSolverType=FinSolverTypes.NELDER_MEAD, tol=1e-8):
        """ Main function to construct the vol surface. """

        s = self._spot_fx_rate
        num_vol_curves = self._num_vol_curves

        if self._volatility_function_type == VolFunctionTypes.CLARK:
            num_parameters = 3
        elif self._volatility_function_type == VolFunctionTypes.SABR:
            num_parameters = 4
        elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_ONE:
            num_parameters = 3
        elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_HALF:
            num_parameters = 3
        elif self._volatility_function_type == VolFunctionTypes.BBG:
            num_parameters = 3
        elif self._volatility_function_type == VolFunctionTypes.CLARK5:
            num_parameters = 5
        else:
            print(self._volatility_function_type)
            raise FinError("Unknown Model Type")

        self._parameters = np.zeros([num_vol_curves, num_parameters])

        num_strikes = 5
        self._strikes = np.zeros([num_vol_curves, num_strikes])
        self._gaps = np.zeros([num_vol_curves, num_strikes])

        self._texp = np.zeros(num_vol_curves)

        self._F0T = np.zeros(num_vol_curves)
        self._rd = np.zeros(num_vol_curves)
        self._rf = np.zeros(num_vol_curves)
        self._K_ATM = np.zeros(num_vol_curves)
        self._deltaATM = np.zeros(num_vol_curves)

        self._K_25D_C = np.zeros(num_vol_curves)
        self._K_25D_P = np.zeros(num_vol_curves)
        self._K_25D_C_MS = np.zeros(num_vol_curves)
        self._K_25D_P_MS = np.zeros(num_vol_curves)
        self._V_25D_MS = np.zeros(num_vol_curves)

        self._K_10D_C = np.zeros(num_vol_curves)
        self._K_10D_P = np.zeros(num_vol_curves)
        self._K_10D_C_MS = np.zeros(num_vol_curves)
        self._K_10D_P_MS = np.zeros(num_vol_curves)
        self._V_10D_MS = np.zeros(num_vol_curves)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        spot_date = self._valuation_date

        for i in range(0, num_vol_curves):

            expiry_date = self._expiry_dates[i]
            texp = (expiry_date - spot_date) / gDaysInYear

            domDF = self._dom_discount_curve._df(texp)
            forDF = self._for_discount_curve._df(texp)
            f = s * forDF/domDF

            self._texp[i] = texp
            self._rd[i] = -np.log(domDF) / texp
            self._rf[i] = -np.log(forDF) / texp
            self._F0T[i] = f

            atm_vol = self._atm_vols[i]

            # This follows exposition in Clarke Page 52
            if self._atmMethod == FinFXATMMethod.SPOT:
                self._K_ATM[i] = s
            elif self._atmMethod == FinFXATMMethod.FWD:
                self._K_ATM[i] = f
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL:
                self._K_ATM[i] = f * np.exp(atm_vol*atm_vol*texp/2.0)
            elif self._atmMethod == FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ:
                self._K_ATM[i] = f * np.exp(-atm_vol*atm_vol*texp/2.0)
            else:
                raise FinError("Unknown Delta Type")

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        ginit = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        x_inits = []
        ginits = []
        for i in range(0, num_vol_curves):

            atm_vol = self._atm_vols[i]

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
            s25 = atm_vol + ms25 + rr25/2.0  # 25D Call
            s50 = atm_vol                   # ATM
            s75 = atm_vol + ms25 - rr25/2.0  # 25D Put (75D Call)

            s10 = atm_vol + ms10 + rr10/2.0  # 10D Call
            s50 = atm_vol                   # ATM
            s90 = atm_vol + ms10 - rr10/2.0  # 10D Put (90D Call)

            if self._volatility_function_type == VolFunctionTypes.CLARK:

                # Our preference is to fit to the 10D wings first
                if self._useMS10DVol is False:
                    # Fit to 25D
                    c0 = np.log(atm_vol)
                    c1 = 2.0 * np.log(s75/s25)
                    c2 = 8.0 * np.log(s25*s75/atm_vol/atm_vol)
                    x_init = [c0, c1, c2]
                else:
                    # Fit to 10D
                    c0 = np.log(atm_vol)
                    c1 = np.log(s90/s10) / 0.80
                    c2 = np.log(s10*s90/atm_vol/atm_vol) / 0.32
                    x_init = [c0, c1, c2]

            elif self._volatility_function_type == VolFunctionTypes.SABR:
                # SABR parameters are alpha, nu, rho
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0
                rho = -0.112
                nu = 0.817

                x_init = [alpha, beta, rho, nu]

            elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_ONE:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0  # FIXED
                rho = -0.112
                nu = 0.817

                x_init = [alpha, rho, nu]

            elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_HALF:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 0.50  # FIXED
                rho = -0.112
                nu = 0.817

                x_init = [alpha, rho, nu]

            elif self._volatility_function_type == VolFunctionTypes.BBG:

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

                x_init = [a, b, c]

            elif self._volatility_function_type == VolFunctionTypes.CLARK5:

                # Our preference is to fit to the 10D wings first
                if self._useMS10DVol is False:
                    # Fit to 25D
                    c0 = np.log(atm_vol)
                    c1 = 2.0 * np.log(s75/s25)
                    c2 = 8.0 * np.log(s25*s75/atm_vol/atm_vol)
                    x_init = [c0, c1, c2, 0.0, 0.0]
                else:
                    # Fit to 10D
                    c0 = np.log(atm_vol)
                    c1 = np.log(s90/s10) / 0.80
                    c2 = np.log(s10*s90/atm_vol/atm_vol) / 0.32
                    x_init = [c0, c1, c2, 0.0, 0.0]

            else:
                raise FinError("Unknown Model Type")

            x_inits.append(x_init)
            ginits.append(ginit)

        delta_method_value = self._deltaMethod.value
        vol_type_value = self._volatility_function_type.value

        for i in range(0, num_vol_curves):

            t = self._texp[i]
            rd = self._rd[i]
            rf = self._rf[i]
            K_ATM = self._K_ATM[i]
            atm_vol = self._atm_vols[i]

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

            res = _solve_to_horizon(s, t, rd, rf,
                                    K_ATM, atm_vol,
                                    ms25DVol, rr25DVol,
                                    ms10DVol, rr10DVol,
                                    delta_method_value, vol_type_value,
                                    self._alpha,
                                    x_inits[i],
                                    ginits[i],
                                    finSolverType,
                                    tol)

            (self._parameters[i, :], self._strikes[i, :], self._gaps[i:],
             self._K_25D_C_MS[i], self._K_25D_P_MS[i],
             self._K_25D_C[i], self._K_25D_P[i],
             self._K_10D_C_MS[i], self._K_10D_P_MS[i],
             self._K_10D_C[i], self._K_10D_P[i]
             ) = res

###############################################################################

    def check_calibration(self, verbose: bool, tol: float = 1e-6):
        """ Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals. """

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valuation_date)
            print("SPOT FX RATE:", self._spot_fx_rate)
            print("ALPHA WEIGHT:", self._alpha)
            print("ATM METHOD:", self._atmMethod)
            print("DELTA METHOD:", self._deltaMethod)
            print("==========================================================")

        K_dummy = 999

        for i in range(0, self._num_vol_curves):

            expiry_date = self._expiry_dates[i]

            if verbose:
                print("TENOR:", self._tenors[i])
                print("EXPIRY DATE:", expiry_date)
                print("IN ATM VOL: %9.6f %%" %
                      (100.0*self._atm_vols[i]))

                if self._useMS25DVol:
                    print("IN MKT STRANGLE 25D VOL: %9.6f %%" %
                          (100.0*self._mktStrangle25DeltaVols[i]))
                    print("IN RSK REVERSAL 25D VOL: %9.6f %%" %
                          (100.0*self._riskReversal25DeltaVols[i]))

                if self._useMS10DVol:
                    print("IN MKT STRANGLE 10D VOL: %9.6f %%" %
                          (100.0*self._mktStrangle10DeltaVols[i]))
                    print("IN RSK REVERSAL 10D VOL: %9.6f %%" %
                          (100.0*self._riskReversal10DeltaVols[i]))

            call = FXVanillaOption(expiry_date,
                                   K_dummy,
                                   self._currency_pair,
                                   OptionTypes.EUROPEAN_CALL,
                                   1.0,
                                   self._notional_currency, )

            put = FXVanillaOption(expiry_date,
                                  K_dummy,
                                  self._currency_pair,
                                  OptionTypes.EUROPEAN_PUT,
                                  1.0,
                                  self._notional_currency)

            ###################################################################
            # AT THE MONEY
            ###################################################################

            if verbose:
                print("==========================================================")
                print("T_(YEARS): ", self._texp[i])
                print("CNT_CPD_RD:%9.6f %%" % (self._rd[i]*100))
                print("CNT_CPD_RF:%9.6f %%" % (self._rf[i]*100))
                print("FWD_RATE:  %9.6f" % (self._F0T[i]))

            sigma_ATM_out = vol_function(self._volatility_function_type.value,
                                         self._parameters[i],
                                         self._strikes[i],
                                         self._gaps[i],
                                         self._F0T[i],
                                         self._K_ATM[i],
                                         self._texp[i])

            if verbose:
                print("==========================================================")
                print("VOL FUNCTION", self._volatility_function_type)
                print("VOL_PARAMETERS:", self._parameters[i])
                print("==========================================================")
                print("OUT_K_ATM:  %9.6f" % (self._K_ATM[i]))
                print("OUT_ATM_VOL: %9.6f %%"
                      % (100.0*sigma_ATM_out))

            diff = sigma_ATM_out - self._atm_vols[i]

            if np.abs(diff) > tol:
                print("FAILED FIT TO ATM VOL IN: %9.6f  OUT: %9.6f  DIFF: %9.6f" %
                      (self._atm_vols[i]*100.0, sigma_ATM_out*100.0,
                       diff * 100.0))

            call._strike_fx_rate = self._K_ATM[i]
            put._strike_fx_rate = self._K_ATM[i]

            model = BlackScholes(sigma_ATM_out)

            delta_call = call.delta(self._valuation_date,
                                    self._spot_fx_rate,
                                    self._dom_discount_curve,
                                    self._for_discount_curve,
                                    model)[self._deltaMethodString]

            delta_put = put.delta(self._valuation_date,
                                  self._spot_fx_rate,
                                  self._dom_discount_curve,
                                  self._for_discount_curve,
                                  model)[self._deltaMethodString]

            if verbose:
                print("CALL_DELTA: % 9.6f  PUT_DELTA: % 9.6f  NET_DELTA: % 9.6f"
                      % (delta_call, delta_put, delta_call + delta_put))

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.25/-0.25
            ###################################################################

            if self._useMS25DVol is True:

                msVol = self._atm_vols[i] + self._mktStrangle25DeltaVols[i]

                if verbose:

                    print("==========================================================")
                    print("MKT STRANGLE 25D VOL IN: %9.6f %%"
                          % (100.0*self._mktStrangle25DeltaVols[i]))

                call._strike_fx_rate = self._K_25D_C_MS[i]
                put._strike_fx_rate = self._K_25D_P_MS[i]

                model = BlackScholes(msVol)

                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
                                      model)[self._deltaMethodString]

                if verbose:
                    print("K_25D_C_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_C_MS[i], 100.0*msVol, delta_call))

                    print("K_25D_P_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_25D_P_MS[i], 100.0*msVol, delta_put))

                call_value = call.value(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)['v']

                put_value = put.value(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                sigma_K_25D_C_MS = vol_function(self._volatility_function_type.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_25D_C_MS[i],
                                                self._texp[i])

                model = BlackScholes(sigma_K_25D_C_MS)
                call_value = call.value(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)['v']

                # THIS IS NOT GOING TO BE 0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                # PUT
                sigma_K_25D_P_MS = vol_function(self._volatility_function_type.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_25D_P_MS[i],
                                                self._texp[i])

                model = BlackScholes(sigma_K_25D_P_MS)
                put_value = put.value(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
                                      model)['v']

                # THIS IS NOT GOING TO BE -0.25 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                    print("FAILED FIT TO 25D MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f" %
                          (mktStrangleValue, mktStrangleValueSkew, diff))

                ###################################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.25, -0.25
                ###################################################################

                call._strike_fx_rate = self._K_25D_C[i]
                put._strike_fx_rate = self._K_25D_P[i]

                sigma_K_25D_C = vol_function(self._volatility_function_type.value,
                                             self._parameters[i],
                                             self._strikes[i],
                                             self._gaps[i],
                                             self._F0T[i],
                                             self._K_25D_C[i],
                                             self._texp[i])

                model = BlackScholes(sigma_K_25D_C)

                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                sigma_K_25D_P = vol_function(self._volatility_function_type.value,
                                             self._parameters[i],
                                             self._strikes[i],
                                             self._gaps[i],
                                             self._F0T[i],
                                             self._K_25D_P[i],
                                             self._texp[i])

                model = BlackScholes(sigma_K_25D_P)

                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                    print("FAILED FIT TO 25D RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f" %
                          (self._riskReversal25DeltaVols[i]*100.0,
                           sigma_RR*100.0,
                           diff*100.0))

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.10/-0.10
            ###################################################################

            if self._useMS10DVol:

                msVol = self._atm_vols[i] + self._mktStrangle10DeltaVols[i]

                if verbose:

                    print("==========================================================")
                    print("MKT STRANGLE 10D VOL IN: %9.6f %%"
                          % (100.0*self._mktStrangle10DeltaVols[i]))

                call._strike_fx_rate = self._K_10D_C_MS[i]
                put._strike_fx_rate = self._K_10D_P_MS[i]

                model = BlackScholes(msVol)

                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
                                      model)[self._deltaMethodString]

                if verbose:
                    print("K_10D_C_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_C_MS[i], 100.0*msVol, delta_call))

                    print("K_10D_P_MS: %9.6f  ATM + MSVOL: %9.6f %%   DELTA: %9.6f"
                          % (self._K_10D_P_MS[i], 100.0*msVol, delta_put))

                call_value = call.value(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)['v']

                put_value = put.value(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                sigma_K_10D_C_MS = vol_function(self._volatility_function_type.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_10D_C_MS[i],
                                                self._texp[i])

                model = BlackScholes(sigma_K_10D_C_MS)
                call_value = call.value(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)['v']

                # THIS IS NOT GOING TO BE 0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                # PUT
                sigma_K_10D_P_MS = vol_function(self._volatility_function_type.value,
                                                self._parameters[i],
                                                self._strikes[i],
                                                self._gaps[i],
                                                self._F0T[i],
                                                self._K_10D_P_MS[i],
                                                self._texp[i])

                model = BlackScholes(sigma_K_10D_P_MS)
                put_value = put.value(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
                                      model)['v']

                # THIS IS NOT GOING TO BE -0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                    print("FAILED FIT TO 10D MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f" %
                          (mktStrangleValue, mktStrangleValueSkew, diff))

                ###################################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.10, -0.10
                ###################################################################

                call._strike_fx_rate = self._K_10D_C[i]
                put._strike_fx_rate = self._K_10D_P[i]

                sigma_K_10D_C = vol_function(self._volatility_function_type.value,
                                             self._parameters[i],
                                             self._strikes[i],
                                             self._gaps[i],
                                             self._F0T[i],
                                             self._K_10D_C[i],
                                             self._texp[i])

                model = BlackScholes(sigma_K_10D_C)

                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(self._valuation_date,
                                        self._spot_fx_rate,
                                        self._dom_discount_curve,
                                        self._for_discount_curve,
                                        model)[self._deltaMethodString]

                sigma_K_10D_P = vol_function(self._volatility_function_type.value,
                                             self._parameters[i],
                                             self._strikes[i],
                                             self._gaps[i],
                                             self._F0T[i],
                                             self._K_10D_P[i],
                                             self._texp[i])

                model = BlackScholes(sigma_K_10D_P)

                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(self._valuation_date,
                                      self._spot_fx_rate,
                                      self._dom_discount_curve,
                                      self._for_discount_curve,
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
                    print("FAILED FIT TO 10D RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f" %
                          (self._riskReversal10DeltaVols[i]*100.0,
                           sigma_RR*100.0,
                           diff*100.0))

###############################################################################

    def implied_dbns(self, lowFX, highFX, numIntervals):
        """ Calculate the pdf for each tenor horizon. Returns a list of
        FinDistribution objects, one for each tenor horizon. """

        dbns = []

        for iTenor in range(0, len(self._tenors)):

            f = self._F0T[iTenor]
            t = self._texp[iTenor]

            dFX = (highFX - lowFX) / numIntervals

            domDF = self._dom_discount_curve._df(t)
            forDF = self._for_discount_curve._df(t)

            rd = -np.log(domDF) / t
            rf = -np.log(forDF) / t

            Ks = []
            vols = []

            for iK in range(0, numIntervals):

                k = lowFX + iK*dFX

                vol = vol_function(self._volatility_function_type.value,
                                   self._parameters[iTenor],
                                   self._strikes[iTenor],
                                   self._gaps[iTenor],
                                   f, k, t)

                Ks.append(k)
                vols.append(vol)

            Ks = np.array(Ks)
            vols = np.array(vols)

            density = option_implied_dbn(
                self._spot_fx_rate, t, rd, rf, Ks, vols)

            dbn = FinDistribution(Ks, density)
            dbns.append(dbn)

        return dbns

###############################################################################

    def plot_vol_curves(self):
        """ Generates a plot of each of the vol discount implied by the market
        and fitted. """

        plt.figure()

        volTypeVal = self._volatility_function_type.value

        for tenorIndex in range(0, self._num_vol_curves):

            atm_vol = self._atm_vols[tenorIndex]*100
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

                sigma = vol_function(volTypeVal, params, strikes, gaps,
                                     f, K, t) * 100.0
                ks.append(K)
                vols.append(sigma)
                K = K + dK

            labelStr = self._tenors[tenorIndex]
            labelStr += " ATM: " + str(atm_vol)[0:6]
            labelStr += " MS25: " + str(msVol25)[0:6]
            labelStr += " RR25: " + str(rrVol25)[0:6]
            labelStr += " MS10: " + str(msVol10)[0:6]
            labelStr += " RR10: " + str(rrVol10)[0:6]

            plt.plot(ks, vols, label=labelStr)
            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = "JNT FIT:" + self._currency_pair + " " +\
                    str(self._volatility_function_type)

            keyStrikes = []
            keyStrikes.append(self._K_ATM[tenorIndex])

            keyVols = []
            for K in keyStrikes:

                sigma = vol_function(volTypeVal, params,
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

                sigma = vol_function(volTypeVal, params,
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
                sigma = vol_function(volTypeVal, params,
                                     strikes, gaps,
                                     f, K, t) * 100.0
                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'ro', markersize=4)

        plt.title(title)
        plt.legend(loc="lower left", bbox_to_anchor=(1, 0))

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUE DATE", self._valuation_date)
        s += label_to_string("FX RATE", self._spot_fx_rate)
        s += label_to_string("CCY PAIR", self._currency_pair)
        s += label_to_string("NOTIONAL CCY", self._notional_currency)
        s += label_to_string("NUM TENORS", self._num_vol_curves)
        s += label_to_string("ATM METHOD", self._atmMethod)
        s += label_to_string("DELTA METHOD", self._deltaMethod)
        s += label_to_string("ALPHA WEIGHT", self._alpha)
        s += label_to_string("VOL FUNCTION", self._volatility_function_type)

        for i in range(0, self._num_vol_curves):

            s += "\n"

            s += label_to_string("TENOR", self._tenors[i])
            s += label_to_string("EXPIRY DATE", self._expiry_dates[i])
            s += label_to_string("TIME (YRS)", self._texp[i])
            s += label_to_string("FWD FX", self._F0T[i])

            s += label_to_string("ATM VOLS", self._atm_vols[i]*100.0)
            s += label_to_string("MS VOLS",
                                 self._mktStrangle25DeltaVols[i]*100.)
            s += label_to_string("RR VOLS",
                                 self._riskReversal25DeltaVols[i]*100.)

            s += label_to_string("ATM Strike", self._K_ATM[i])
            s += label_to_string("ATM Delta", self._deltaATM[i])

            s += label_to_string("K_ATM", self._K_ATM[i])

            s += label_to_string("MS 25D Call Strike", self._K_25D_C_MS[i])
            s += label_to_string("MS 25D Put Strike", self._K_25D_P_MS[i])
            s += label_to_string("SKEW 25D CALL STRIKE", self._K_25D_C[i])
            s += label_to_string("SKEW 25D PUT STRIKE", self._K_25D_P[i])
            s += label_to_string("PARAMS", self._parameters[i])

            s += label_to_string("MS 10D Call Strike", self._K_10D_C_MS[i])
            s += label_to_string("MS 10D Put Strike", self._K_10D_P_MS[i])
            s += label_to_string("SKEW 10D CALL STRIKE", self._K_10D_C[i])
            s += label_to_string("SKEW 10D PUT STRIKE", self._K_10D_P[i])

        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
