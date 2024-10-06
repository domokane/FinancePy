##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt
from numba import njit, float64, int64

from ...utils.error import FinError
from ...utils.date import Date
from ...utils.global_vars import g_days_in_year
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
from ...models.volatility_fns import VolFuncTypes
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
    """This is the objective function used in the determination of the FX
    option implied strike which is computed in the class below."""

    s = args[0]
    t = args[1]
    r_d = args[2]
    r_f = args[3]
    volatility = args[4]
    delta_method_value = args[5]
    option_type_value = args[6]
    delta_target = args[7]

    delta_out = fast_delta(
        s, t, K, r_d, r_f, volatility, delta_method_value, option_type_value
    )

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
        if k > strikes[i - 1] and k <= strikes[i]:
            index = i
            break

    if index == 0:
        raise FinError("Value not bracketed")

    k0 = strikes[index - 1]
    k1 = strikes[index]
    v0 = gaps[index - 1]
    v1 = gaps[index]
    v = ((k - k0) * v1 + (k1 - k) * v0) / (k1 - k0)
    return v


###############################################################################
# Do not cache this function


@njit(fastmath=True)  # , cache=True)
def _obj(params, *args):
    """Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the vol_type_value
    """

    s = args[0]
    t = args[1]
    r_d = args[2]
    r_f = args[3]
    k_atm = args[4]
    atm_vol = args[5]

    k_25d_c_ms = args[6]
    k_25d_p_ms = args[7]
    v_25d_ms_target = args[8]
    target_25d_rr_vol = args[9]

    k_10d_c_ms = args[10]
    k_10d_p_ms = args[11]
    v_10d_ms_target = args[12]
    target_10d_rr_vol = args[13]

    delta_method_value = args[14]
    vol_type_value = args[15]
    alpha = args[16]

    strikes_null = np.zeros(1)
    gaps_null = np.zeros(1)

    f = s * np.exp((r_d - r_f) * t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atm_curve_vol = vol_function(
        vol_type_value, params, strikes_null, gaps_null, f, k_atm, t
    )

    term_atm = (atm_vol - atm_curve_vol) ** 2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25d strikes
    ###########################################################################

    if target_25d_rr_vol > -999.0:

        sigma_k_25d_c_ms = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_25d_c_ms, t
        )

        v_25d_c_ms = bs_value(
            s,
            t,
            k_25d_c_ms,
            r_d,
            r_f,
            sigma_k_25d_c_ms,
            OptionTypes.EUROPEAN_CALL.value,
        )

        sigma_k_25d_p_ms = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_25d_p_ms, t
        )

        v_25d_p_ms = bs_value(
            s,
            t,
            k_25d_p_ms,
            r_d,
            r_f,
            sigma_k_25d_p_ms,
            OptionTypes.EUROPEAN_PUT.value,
        )

        v_25d_ms = v_25d_c_ms + v_25d_p_ms
        term_25d_1 = (v_25d_ms - v_25d_ms_target) ** 2

    else:

        term_25d_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target_25d_rr_vol > -999.0:

        k_25d_c = _solver_for_smile_strike(
            s,
            t,
            r_d,
            r_f,
            OptionTypes.EUROPEAN_CALL.value,
            vol_type_value,
            +0.2500,
            delta_method_value,
            k_25d_c_ms,
            params,
            strikes_null,
            gaps_null,
        )

        sigma_k_25d_c = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_25d_c, t
        )

        k_25d_p = _solver_for_smile_strike(
            s,
            t,
            r_d,
            r_f,
            OptionTypes.EUROPEAN_PUT.value,
            vol_type_value,
            -0.2500,
            delta_method_value,
            k_25d_p_ms,
            params,
            strikes_null,
            gaps_null,
        )

        sigma_k_25d_p = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_25d_p, t
        )

        sigma_25d_rr = sigma_k_25d_c - sigma_k_25d_p
        term_25d_2 = (sigma_25d_rr - target_25d_rr_vol) ** 2

    else:

        term_25d_2 = 0.0

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 10d strikes
    ###########################################################################

    if target_10d_rr_vol > -999.0:

        sigma_k_10d_c_ms = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_10d_c_ms, t
        )

        v_10d_c_ms = bs_value(
            s,
            t,
            k_10d_c_ms,
            r_d,
            r_f,
            sigma_k_10d_c_ms,
            OptionTypes.EUROPEAN_CALL.value,
        )

        sigma_k_10d_p_ms = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_10d_p_ms, t
        )

        v_10d_p_ms = bs_value(
            s,
            t,
            k_10d_p_ms,
            r_d,
            r_f,
            sigma_k_10d_p_ms,
            OptionTypes.EUROPEAN_PUT.value,
        )

        v_10d_ms = v_10d_c_ms + v_10d_p_ms
        term10d_1 = (v_10d_ms - v_10d_ms_target) ** 2

    else:

        term10d_1 = 0.0

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    if target_10d_rr_vol > -999.0:

        k_10d_c = _solver_for_smile_strike(
            s,
            t,
            r_d,
            r_f,
            OptionTypes.EUROPEAN_CALL.value,
            vol_type_value,
            +0.1000,
            delta_method_value,
            k_10d_c_ms,
            params,
            strikes_null,
            gaps_null,
        )

        sigma_k_10d_c = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_10d_c, t
        )

        k_10d_p = _solver_for_smile_strike(
            s,
            t,
            r_d,
            r_f,
            OptionTypes.EUROPEAN_PUT.value,
            vol_type_value,
            -0.1000,
            delta_method_value,
            k_10d_p_ms,
            params,
            strikes_null,
            gaps_null,
        )

        sigma_k_10d_p = vol_function(
            vol_type_value, params, strikes_null, gaps_null, f, k_10d_p, t
        )

        sigma_10d_rr = sigma_k_10d_c - sigma_k_10d_p
        term10d_2 = (sigma_10d_rr - target_10d_rr_vol) ** 2

    else:

        term10d_2 = 0.0

    ###########################################################################
    # Alpha interpolates between fitting only ATM and 25d when alpha = 0.0 and
    # fitting only ATM and 10d when alpha = 1.0. Equal when alpha = 0.50.
    ###########################################################################

    tot = term_atm
    tot = tot + (1.0 - alpha) * (term_25d_1 + term_25d_2)
    tot = tot + alpha * (term10d_1 + term10d_2)
    return tot


###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################

# THIS FUNCTION IS NOT USED CURRENTLY


# @njit(fastmath=True)  # , cache=True)
def _obj_gap(gaps, *args):
    """Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the vol_type_value
    """

    s = args[0]
    t = args[1]
    r_d = args[2]
    r_f = args[3]
    k_atm = args[4]
    atm_vol = args[5]

    k_25d_c_ms = args[6]
    k_25d_p_ms = args[7]
    v_25d_ms_target = args[8]
    target_25d_rr_vol = args[9]

    k_10d_c_ms = args[10]
    k_10d_p_ms = args[11]
    v_10d_ms_target = args[12]
    target_10d_rr_vol = args[13]

    delta_method_value = args[14]
    vol_type_value = args[15]
    params = args[16]

    strikes = [k_10d_p_ms, k_25d_p_ms, k_atm, k_25d_c_ms, k_10d_c_ms]
    strikes = np.array(strikes)

    f = s * np.exp((r_d - r_f) * t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atm_curve_vol = vol_function(
        vol_type_value, params, strikes, gaps, f, k_atm, t
    )

    print("atm_curve_vol", atm_curve_vol)

    term_atm = (atm_vol - atm_curve_vol) ** 2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 25d strikes
    ###########################################################################

    sigma_k_25d_c_ms = vol_function(
        vol_type_value, params, strikes, gaps, f, k_25d_c_ms, t
    )

    print("sigma_k_25d_c_ms", sigma_k_25d_c_ms)

    v_25d_c_ms = bs_value(
        s,
        t,
        k_25d_c_ms,
        r_d,
        r_f,
        sigma_k_25d_c_ms,
        OptionTypes.EUROPEAN_CALL.value,
    )

    sigma_k_25d_p_ms = vol_function(
        vol_type_value, params, strikes, gaps, f, k_25d_p_ms, t
    )

    print("sigma_k_25d_p_ms", sigma_k_25d_p_ms)

    v_25d_p_ms = bs_value(
        s,
        t,
        k_25d_p_ms,
        r_d,
        r_f,
        sigma_k_25d_p_ms,
        OptionTypes.EUROPEAN_PUT.value,
    )

    v_25d_ms = v_25d_c_ms + v_25d_p_ms
    term_25d_1 = (v_25d_ms - v_25d_ms_target) ** 2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    k_25d_c = _solver_for_smile_strike(
        s,
        t,
        r_d,
        r_f,
        OptionTypes.EUROPEAN_CALL.value,
        vol_type_value,
        +0.2500,
        delta_method_value,
        k_25d_c_ms,
        params,
        strikes,
        gaps,
    )

    sigma_k_25d_c = vol_function(
        vol_type_value, params, strikes, gaps, f, k_25d_c, t
    )

    print("sigma_k_25d_c", sigma_k_25d_c)

    k_25d_p = _solver_for_smile_strike(
        s,
        t,
        r_d,
        r_f,
        OptionTypes.EUROPEAN_PUT.value,
        vol_type_value,
        -0.2500,
        delta_method_value,
        k_25d_p_ms,
        params,
        strikes,
        gaps,
    )

    sigma_k_25d_p = vol_function(
        vol_type_value, params, strikes, gaps, f, k_25d_p, t
    )

    print("sigma_k_25d_p", sigma_k_25d_p)

    sigma_25d_rr = sigma_k_25d_c - sigma_k_25d_p
    term_25d_2 = (sigma_25d_rr - target_25d_rr_vol) ** 2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS 10d strikes
    ###########################################################################

    sigma_k_10d_c_ms = vol_function(
        vol_type_value, params, strikes, gaps, f, k_10d_c_ms, t
    )

    print("sigma_k_10d_c_ms", sigma_k_10d_c_ms)

    v_10d_c_ms = bs_value(
        s,
        t,
        k_10d_c_ms,
        r_d,
        r_f,
        sigma_k_10d_c_ms,
        OptionTypes.EUROPEAN_CALL.value,
    )

    sigma_k_10d_p_ms = vol_function(
        vol_type_value, params, strikes, gaps, f, k_10d_p_ms, t
    )

    print("sigma_k_10d_p_ms", sigma_k_10d_p_ms)

    v_10d_p_ms = bs_value(
        s,
        t,
        k_10d_p_ms,
        r_d,
        r_f,
        sigma_k_10d_p_ms,
        OptionTypes.EUROPEAN_PUT.value,
    )

    v_10d_ms = v_10d_c_ms + v_10d_p_ms
    term10d_1 = (v_10d_ms - v_10d_ms_target) ** 2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    k_10d_c = _solver_for_smile_strike(
        s,
        t,
        r_d,
        r_f,
        OptionTypes.EUROPEAN_CALL.value,
        vol_type_value,
        +0.1000,
        delta_method_value,
        k_10d_c_ms,
        params,
        strikes,
        gaps,
    )

    sigma_k_10d_c = vol_function(
        vol_type_value, params, strikes, gaps, f, k_10d_c, t
    )

    print("SIGMA_k_10d_c", sigma_k_10d_c)

    print("INIT k_10d_p_ms", k_10d_p_ms)

    k_10d_p = _solver_for_smile_strike(
        s,
        t,
        r_d,
        r_f,
        OptionTypes.EUROPEAN_PUT.value,
        vol_type_value,
        -0.1000,
        delta_method_value,
        k_10d_p_ms,
        params,
        strikes,
        gaps,
    )

    print("k_10d_p", k_10d_p)
    sigma_k_10d_p = vol_function(
        vol_type_value, params, strikes, gaps, f, k_10d_p, t
    )

    print("SIGMA_k_10d_p", sigma_k_10d_p)

    sigma_10d_rr = sigma_k_10d_c - sigma_k_10d_p
    term10d_2 = (sigma_10d_rr - target_10d_rr_vol) ** 2

    ###########################################################################
    # Alpha interpolates between fitting only ATM and 25d when alpha = 0.0 and
    # fitting only ATM and 10d when alpha = 1.0. Equal when alpha = 0.50.
    ###########################################################################

    tot = term_atm
    tot = tot + (term_25d_1 + term_25d_2)
    tot = tot + (term10d_1 + term10d_2)
    return tot


###############################################################################


def _solve_to_horizon(
    s,
    t,
    rd,
    rf,
    k_atm,
    atm_vol,
    ms_25d_vol,
    rr_25d_vol,
    ms_10d_vol,
    rr_10d_vol,
    delta_method_value,
    vol_type_value,
    alpha,
    x_inits,
    ginits,
    fin_solver_type,
    tol,
):

    ###########################################################################
    # Determine the price of a market strangle from market strangle
    # Need to price a call and put that agree with market strangle
    ###########################################################################

    use10d = True
    use25d = True

    if ms_25d_vol == -999.0:
        use25d = False

    if ms_10d_vol == -999.0:
        use10d = False

    if use25d is True:

        vol_25d_ms = atm_vol + ms_25d_vol

        k_25d_c_ms = solve_for_strike(
            s,
            t,
            rd,
            rf,
            OptionTypes.EUROPEAN_CALL.value,
            +0.2500,
            delta_method_value,
            vol_25d_ms,
        )

        k_25d_p_ms = solve_for_strike(
            s,
            t,
            rd,
            rf,
            OptionTypes.EUROPEAN_PUT.value,
            -0.2500,
            delta_method_value,
            vol_25d_ms,
        )

        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        v_25d_c_ms = bs_value(
            s,
            t,
            k_25d_c_ms,
            rd,
            rf,
            vol_25d_ms,
            OptionTypes.EUROPEAN_CALL.value,
        )

        v_25d_p_ms = bs_value(
            s,
            t,
            k_25d_p_ms,
            rd,
            rf,
            vol_25d_ms,
            OptionTypes.EUROPEAN_PUT.value,
        )

        # Market price of strangle in the domestic currency
        v_25d_ms = v_25d_c_ms + v_25d_p_ms

    else:

        vol_25d_ms = -999.0
        k_25d_c_ms = 0.0
        k_25d_p_ms = 0.0
        v_25d_c_ms = 0.0
        v_25d_p_ms = 0.0
        v_25d_ms = 0.0

    ###########################################################################

    if use10d is True:

        vol_10d_ms = atm_vol + ms_10d_vol

        k_10d_c_ms = solve_for_strike(
            s,
            t,
            rd,
            rf,
            OptionTypes.EUROPEAN_CALL.value,
            +0.1000,
            delta_method_value,
            vol_10d_ms,
        )

        k_10d_p_ms = solve_for_strike(
            s,
            t,
            rd,
            rf,
            OptionTypes.EUROPEAN_PUT.value,
            -0.1000,
            delta_method_value,
            vol_10d_ms,
        )

        # USE MARKET STRANGLE VOL TO DETERMINE PRICE OF A MARKET STRANGLE
        v_10d_c_ms = bs_value(
            s,
            t,
            k_10d_c_ms,
            rd,
            rf,
            vol_10d_ms,
            OptionTypes.EUROPEAN_CALL.value,
        )

        v_10d_p_ms = bs_value(
            s,
            t,
            k_10d_p_ms,
            rd,
            rf,
            vol_10d_ms,
            OptionTypes.EUROPEAN_PUT.value,
        )

        # Market price of strangle in the domestic currency
        v_10d_ms = v_10d_c_ms + v_10d_p_ms

    else:

        vol_10d_ms = -999.0
        k_10d_c_ms = 0.0
        k_10d_p_ms = 0.0
        v_10d_c_ms = 0.0
        v_10d_p_ms = 0.0
        v_10d_ms = 0.0

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    # tol = 1e-8

    args = (
        s,
        t,
        rd,
        rf,
        k_atm,
        atm_vol,
        k_25d_c_ms,
        k_25d_p_ms,
        v_25d_ms,
        rr_25d_vol,
        k_10d_c_ms,
        k_10d_p_ms,
        v_10d_ms,
        rr_10d_vol,
        delta_method_value,
        vol_type_value,
        alpha,
    )

    # Nelder-Mead (both SciPy & Numba) is quicker, but occasionally fails
    # to converge, so for those cases try again with CG
    # Numba version is quicker, but can be slightly away from CG output
    try:
        if fin_solver_type == FinSolverTypes.NELDER_MEAD_NUMBA:
            xopt = nelder_mead(
                _obj,
                np.array(x_inits),
                bounds=np.array([[], []]).T,
                args=args,
                tol_f=tol,
                tol_x=tol,
                max_iter=1000,
            )
        elif fin_solver_type == FinSolverTypes.NELDER_MEAD:
            opt = minimize(_obj, x_inits, args, method="Nelder-Mead", tol=tol)
            xopt = opt.x
        elif fin_solver_type == FinSolverTypes.CONJUGATE_GRADIENT:
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)
            xopt = opt.x
    except Exception:
        # If convergence fails try again with CG if necessary
        if fin_solver_type != FinSolverTypes.CONJUGATE_GRADIENT:
            print("Failed to converge, will try CG")
            opt = minimize(_obj, x_inits, args, method="CG", tol=tol)

            xopt = opt.x

    params = np.array(xopt)

    strikes = [k_10d_p_ms, k_25d_p_ms, k_atm, k_10d_c_ms, k_25d_c_ms]
    strikes = np.array(strikes)
    gaps = np.zeros(5)

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    if 1 == 0:

        tol = 1e-12

        args = (
            s,
            t,
            rd,
            rf,
            k_atm,
            atm_vol,
            k_25d_c_ms,
            k_25d_p_ms,
            v_25d_ms,
            rr_25d_vol,
            k_10d_c_ms,
            k_10d_p_ms,
            v_10d_ms,
            rr_10d_vol,
            delta_method_value,
            vol_type_value,
            params,
        )

        opt = minimize(_obj_gap, ginits, args, method="Nelder-Mead", tol=tol)
        xopt = opt.x
        gaps = np.array(xopt)

        print("SOLVED")

    # Removed this as it causes discontinuity
    #    f = s * np.exp((rd-rf)*t)
    #    interpATMVol = vol_function(vol_type_value, params,
    #                                   strikes, gaps, f, k_atm, t)

    #    diff = atm_vol - interpATMVol
    #    gaps[2] = diff

    ###########################################################################

    if use25d is False:
        k_25d_c_ms = k_atm
        k_25d_p_ms = k_atm

    k_25d_c = _solver_for_smile_strike(
        s,
        t,
        rd,
        rf,
        OptionTypes.EUROPEAN_CALL.value,
        vol_type_value,
        +0.2500,
        delta_method_value,
        k_25d_c_ms,
        params,
        strikes,
        gaps,
    )

    k_25d_p = _solver_for_smile_strike(
        s,
        t,
        rd,
        rf,
        OptionTypes.EUROPEAN_PUT.value,
        vol_type_value,
        -0.2500,
        delta_method_value,
        k_25d_p_ms,
        params,
        strikes,
        gaps,
    )

    if use10d is False:
        k_10d_c_ms = k_atm
        k_10d_p_ms = k_atm

    k_10d_c = _solver_for_smile_strike(
        s,
        t,
        rd,
        rf,
        OptionTypes.EUROPEAN_CALL.value,
        vol_type_value,
        +0.1000,
        delta_method_value,
        k_10d_c_ms,
        params,
        strikes,
        gaps,
    )

    k_10d_p = _solver_for_smile_strike(
        s,
        t,
        rd,
        rf,
        OptionTypes.EUROPEAN_PUT.value,
        vol_type_value,
        -0.1000,
        delta_method_value,
        k_10d_p_ms,
        params,
        strikes,
        gaps,
    )

    return (
        params,
        strikes,
        gaps,
        k_25d_c_ms,
        k_25d_p_ms,
        k_25d_c,
        k_25d_p,
        k_10d_c_ms,
        k_10d_p_ms,
        k_10d_c,
        k_10d_p,
    )


###############################################################################


@njit(
    float64(
        int64, float64[:], float64[:], float64[:], float64, float64, float64
    ),
    cache=True,
    fastmath=True,
)
def vol_function(vol_function_type_value, params, strikes, gaps, f, k, t):
    """Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book."""

    #    print("vol_function", vol_function_type_value)

    if len(strikes) == 1:
        gap_k = 0.0
    else:
        gap_k = _interpolate_gap(k, strikes, gaps)

    if vol_function_type_value == VolFuncTypes.CLARK.value:
        vol = vol_function_clark(params, f, k, t) + gap_k
        return vol

    if vol_function_type_value == VolFuncTypes.SABR.value:
        vol = vol_function_sabr(params, f, k, t) + gap_k
        return vol

    if vol_function_type_value == VolFuncTypes.SABR_BETA_HALF.value:
        vol = vol_function_sabr_beta_half(params, f, k, t) + gap_k
        return vol

    if vol_function_type_value == VolFuncTypes.SABR_BETA_ONE.value:
        vol = vol_function_sabr_beta_one(params, f, k, t) + gap_k
        return vol

    if vol_function_type_value == VolFuncTypes.BBG.value:
        vol = vol_function_bloomberg(params, f, k, t) + gap_k
        return vol

    if vol_function_type_value == VolFuncTypes.CLARK5.value:
        vol = vol_function_clark(params, f, k, t) + gap_k
        return vol

    raise FinError("Unknown Model Type")


###############################################################################


@njit(cache=True, fastmath=True)
def _delta_fit(k, *args):
    """This is the objective function used in the determination of the FX
    Option implied strike which is computed in the class below. I map it into
    inverse normcdf space to avoid the flat slope of this function at low vol
    and high K. It speeds up the code as it allows initial values close to
    the solution to be used."""

    vol_type_value = args[0]
    s = args[1]
    t = args[2]
    r_d = args[3]
    r_f = args[4]
    option_type_value = args[5]
    delta_type_value = args[6]
    inverse_delta_target = args[7]
    params = args[8]
    strikes = args[9]
    gaps = args[10]

    f = s * np.exp((r_d - r_f) * t)
    v = vol_function(vol_type_value, params, strikes, gaps, f, k, t)

    delta_out = fast_delta(
        s, t, k, r_d, r_f, v, delta_type_value, option_type_value
    )

    inverse_delta_out = norminvcdf(np.abs(delta_out))
    inv_obj_fn = inverse_delta_target - inverse_delta_out

    #    print(k, f, v, delta_out, inv_obj_fn)

    return inv_obj_fn


###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        int64,
        int64,
        float64,
        int64,
        float64,
        float64[:],
        float64[:],
        float64[:],
    ),
    fastmath=True,
    cache=False,
)
def _solver_for_smile_strike(
    s,
    t,
    rd,
    rf,
    option_type_value,
    vol_type_value,
    delta_target,
    delta_method_value,
    initial_guess,
    parameters,
    strikes,
    gaps,
):
    """Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike."""

    inverse_delta_target = norminvcdf(np.abs(delta_target))

    argtuple = (
        vol_type_value,
        s,
        t,
        rd,
        rf,
        option_type_value,
        delta_method_value,
        inverse_delta_target,
        parameters,
        strikes,
        gaps,
    )

    K = newton_secant(
        _delta_fit, x0=initial_guess, args=argtuple, tol=1e-8, maxiter=50
    )

    return K


###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


@njit(
    float64(
        float64, float64, float64, float64, int64, float64, int64, float64
    ),
    fastmath=True,
)
def solve_for_strike(
    spot_fx_rate,
    t_del,
    rd,
    rf,
    option_type_value,
    delta_target,
    delta_method_value,
    volatility,
):
    """This function determines the implied strike of an FX option
    given a delta and the other option details. It uses a one-dimensional
    Newton root search algorithm to determine the strike that matches an
    input volatility."""

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

        dom_df = np.exp(-rd * t_del)
        for_df = np.exp(-rf * t_del)

        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        fwd_fx_rate = spot_fx_rate * for_df / dom_df
        vsqrtt = volatility * np.sqrt(t_del)
        arg = delta_target * phi / for_df  # CHECK THIS !!!
        norminvdelta = norminvcdf(arg)
        K = fwd_fx_rate * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt / 2.0))
        return K

    if delta_method_value == FinFXDeltaMethod.FORWARD_DELTA.value:

        dom_df = np.exp(-rd * t_del)
        for_df = np.exp(-rf * t_del)

        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            phi = +1.0
        else:
            phi = -1.0

        fwd_fx_rate = spot_fx_rate * for_df / dom_df
        vsqrtt = volatility * np.sqrt(t_del)
        arg = delta_target * phi
        norminvdelta = norminvcdf(arg)
        K = fwd_fx_rate * np.exp(-vsqrtt * (phi * norminvdelta - vsqrtt / 2.0))
        return K

    if delta_method_value == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

        argtuple = (
            spot_fx_rate,
            t_del,
            rd,
            rf,
            volatility,
            delta_method_value,
            option_type_value,
            delta_target,
        )

        K = newton_secant(
            _g, x0=spot_fx_rate, args=argtuple, tol=1e-7, maxiter=50
        )

        return K

    if delta_method_value == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

        argtuple = (
            spot_fx_rate,
            t_del,
            rd,
            rf,
            volatility,
            delta_method_value,
            option_type_value,
            delta_target,
        )

        K = newton_secant(
            _g, x0=spot_fx_rate, args=argtuple, tol=1e-7, maxiter=50
        )

        return K

    raise FinError("Unknown FinFXDeltaMethod")


###############################################################################


class FXVolSurfacePlus:
    """Class to perform a calibration of a chosen parametrised surface to the
    prices of FX options at different strikes and expiry tenors. The
    calibration inputs are the ATM and 25 and 10 Delta volatilities in terms of
    the market strangle amd risk reversals. There is a choice of volatility
    function from cubic in delta to full SABR. Check out VolFuncTypes.
    Parameter alpha [0,1] is used to interpolate between fitting only 25d when
    alpha=0 to fitting only 10d when alpha=1.0. Alpha=0.5 assigns equal weights
    A vol function with more parameters will give a better fit. Of course. But
    it might also overfit. Visualising the volatility curve is useful. Also,
    there is no guarantee that the implied pdf will be positive."""

    def __init__(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        currency_pair: str,
        notional_currency: str,
        domestic_curve: DiscountCurve,
        foreign_curve: DiscountCurve,
        tenors: list,
        atm_vols: (list, np.ndarray),
        ms_25_delta_vols: (list, np.ndarray),
        rr_25_delta_vols: (list, np.ndarray),
        ms_10_delta_vols: (list, np.ndarray),
        rr_10_delta_vols: (list, np.ndarray),
        alpha: float,
        atm_method: FinFXATMMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL,
        delta_method: FinFXDeltaMethod = FinFXDeltaMethod.SPOT_DELTA,
        vol_func_type: VolFuncTypes = VolFuncTypes.CLARK,
        fin_solver_type: FinSolverTypes = FinSolverTypes.NELDER_MEAD,
        tol: float = 1e-8,
    ):
        """Create the FinFXVolSurfacePlus object by passing in market vol data
        for ATM, 25 Delta and 10 Delta strikes. The alpha weight shifts the
        fitting between 25d and 10d. Alpha = 0.0 is 100% 25d while alpha = 1.0
        is 100% 10d. An alpha of 0.50 is equally weighted."""

        # I want to allow Nones for some of the market inputs
        if ms_10_delta_vols is None:
            ms_10_delta_vols = []

        if rr_10_delta_vols is None:
            rr_10_delta_vols = []

        if ms_25_delta_vols is None:
            ms_25_delta_vols = []

        if rr_25_delta_vols is None:
            rr_25_delta_vols = []

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt
        self.spot_fx_rate = spot_fx_rate
        self.currency_pair = currency_pair

        if len(currency_pair) != 6:
            raise FinError("Currency pair must be 6 characters.")

        self.for_name = self.currency_pair[0:3]
        self.dom_name = self.currency_pair[3:6]

        self.notional_currency = notional_currency
        self.domestic_curve = domestic_curve
        self.foreign_curve = foreign_curve
        self.num_vol_curves = len(tenors)
        self.tenors = tenors

        if len(atm_vols) != self.num_vol_curves:
            raise FinError("Number ATM vols must equal number of tenors")

        self.atm_vols = np.array(atm_vols) / 100.0

        self.use_ms_25d_vol = True
        self.use_rr_25d_vol = True
        self.use_ms_10d_vol = True
        self.use_rr_10d_vol = True

        # Some of these can be missing which is signified by length zero
        n = len(ms_25_delta_vols)

        if n != self.num_vol_curves and n != 0:
            raise FinError("Number MS25d vols must equal number of tenors")

        if n == 0:
            self.use_ms_25d_vol = False

        n = len(rr_25_delta_vols)

        if n != self.num_vol_curves and n != 0:
            raise FinError("Number RR25d vols must equal number of tenors")

        if n == 0:
            self.use_rr_25d_vol = False

        n = len(ms_10_delta_vols)

        if n != self.num_vol_curves and n != 0:
            raise FinError("Number MS10d vols must equal number of tenors")

        if n == 0:
            self.use_ms_10d_vol = False

        n = len(rr_10_delta_vols)

        if n != self.num_vol_curves and n != 0:
            raise FinError("Number RR10d vols must equal number of tenors")

        if n == 0:
            self.use_rr_10d_vol = False

        if self.use_ms_10d_vol != self.use_rr_10d_vol:
            raise FinError("You must provide both 10d RR + 10d MS or neither")

        if self.use_ms_25d_vol != self.use_rr_25d_vol:
            raise FinError("You must provide both 25d RR + 25d MS or neither")

        if self.use_ms_10d_vol is False and self.use_ms_25d_vol is False:
            raise FinError(
                "No MS and RR. You must provide 10d or 25d MS + RR."
            )

        self.ms_25_delta_vols = np.array(ms_25_delta_vols) / 100.0
        self.rr_25_delta_vols = np.array(rr_25_delta_vols) / 100.0
        self.ms_10_delta_vols = np.array(ms_10_delta_vols) / 100.0
        self.rr_10_delta_vols = np.array(rr_10_delta_vols) / 100.0

        if alpha < 0.0 or alpha > 1.0:
            raise FinError("Alpha must be between 0.0 and 1.0")

        self.alpha = alpha

        self.atm_method = atm_method
        self.delta_method = delta_method

        if self.delta_method == FinFXDeltaMethod.SPOT_DELTA:
            self.delta_method_string = "pips_spot_delta"
        elif self.delta_method == FinFXDeltaMethod.FORWARD_DELTA:
            self.delta_method_string = "pips_fwd_delta"
        elif self.delta_method == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ:
            self.delta_method_string = "pct_spot_delta_prem_adj"
        elif self.delta_method == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ:
            self.delta_method_string = "pct_fwd_delta_prem_adj"
        else:
            raise FinError("Unknown Delta Type")

        self.vol_func_type = vol_func_type
        self.tenor_index = 0

        self.expiry_dts = []
        for i in range(0, self.num_vol_curves):
            expiry_dt = value_dt.add_tenor(tenors[i])
            self.expiry_dts.append(expiry_dt)

        self._build_vol_surface(fin_solver_type=fin_solver_type, tol=tol)

    ###########################################################################

    def vol_from_strike_date(self, K, expiry_dt):
        """Interpolates the Black-Scholes volatility from the volatility
        surface given call option strike and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be
        overriden by a provided delta convention. The resulting volatilities
        are then determined for each bracketing expiry time and linear
        interpolation is done in variance space and then converted back to a
        lognormal volatility."""

        t_exp = (expiry_dt - self.value_dt) / g_days_in_year

        vol_type_value = self.vol_func_type.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self.num_vol_curves

        if num_curves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif t_exp <= self.t_exp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif t_exp >= self.t_exp[-1]:

            index0 = len(self.t_exp) - 1
            index1 = len(self.t_exp) - 1

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if t_exp <= self.t_exp[i] and t_exp > self.t_exp[i - 1]:
                    index0 = i - 1
                    index1 = i
                    break

        fwd0 = self.fwd[index0]
        fwd1 = self.fwd[index1]

        t0 = self.t_exp[index0]
        t1 = self.t_exp[index1]

        vol0 = vol_function(
            vol_type_value,
            self.parameters[index0],
            self.strikes[index0],
            self.gaps[index0],
            fwd0,
            K,
            t0,
        )

        if index1 != index0:

            vol1 = vol_function(
                vol_type_value,
                self.parameters[index1],
                self.strikes[index1],
                self.gaps[index1],
                fwd1,
                K,
                t1,
            )

        else:

            vol1 = vol0

        # In the expiry time dimension, both volatilities are interpolated
        # at the same strikes but different deltas.
        vart0 = vol0 * vol0 * t0
        vart1 = vol1 * vol1 * t1

        if np.abs(t1 - t0) > 1e-6:
            vart = ((t_exp - t0) * vart1 + (t1 - t_exp) * vart0) / (t1 - t0)

            if vart < 0.0:
                raise FinError("Negative variance.")

            volt = np.sqrt(vart / t_exp)

        else:
            volt = vol1

        return volt

    ###########################################################################

    def delta_to_strike(self, call_delta, expiry_dt, delta_method):
        """Interpolates the strike at a delta and expiry date. Linear
        time to expiry interpolation is used in strike."""

        t_exp = (expiry_dt - self.value_dt) / g_days_in_year

        vol_type_value = self.vol_func_type.value

        s = self.spot_fx_rate

        if delta_method is None:
            delta_method_value = self.delta_method.value
        else:
            delta_method_value = delta_method.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self.num_vol_curves

        # If there is only one time horizon then assume flat vol to this time
        if num_curves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif t_exp <= self.t_exp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif t_exp > self.t_exp[-1]:

            index0 = len(self.t_exp) - 1
            index1 = len(self.t_exp) - 1

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if t_exp <= self.t_exp[i] and t_exp > self.t_exp[i - 1]:
                    index0 = i - 1
                    index1 = i
                    break

        #######################################################################

        t0 = self.t_exp[index0]
        t1 = self.t_exp[index1]

        initial_guess = self.k_atm[index0]

        k0 = _solver_for_smile_strike(
            s,
            t_exp,
            self.rd[index0],
            self.rf[index0],
            OptionTypes.EUROPEAN_CALL.value,
            vol_type_value,
            call_delta,
            delta_method_value,
            initial_guess,
            self.parameters[index0],
            self.strikes[index0],
            self.gaps[index0],
        )

        if index1 != index0:

            k1 = _solver_for_smile_strike(
                s,
                t_exp,
                self.rd[index1],
                self.rf[index1],
                OptionTypes.EUROPEAN_CALL.value,
                vol_type_value,
                call_delta,
                delta_method_value,
                initial_guess,
                self.parameters[index1],
                self.strikes[index1],
                self.gaps[index1],
            )
        else:

            k1 = k0

        # In the expiry time dimension, both volatilities are interpolated
        # at the same strikes but different deltas.

        if np.abs(t1 - t0) > 1e-6:

            K = ((t_exp - t0) * k1 + (t1 - t_exp) * k0) / (t1 - t0)

        else:

            K = k1

        return K

    ###########################################################################

    def vol_from_delta_date(self, call_delta, expiry_dt, delta_method=None):
        """Interpolates the Black-Scholes volatility from the volatility
        surface given a call option delta and expiry date. Linear interpolation
        is done in variance space. The smile strikes at bracketed dates are
        determined by determining the strike that reproduces the provided delta
        value. This uses the calibration delta convention, but it can be
        overriden by a provided delta convention. The resulting volatilities
        are then determined for each bracketing expiry time and linear
        interpolation is done in variance space and then converted back to a
        lognormal volatility."""

        t_exp = (expiry_dt - self.value_dt) / g_days_in_year

        vol_type_value = self.vol_func_type.value

        s = self.spot_fx_rate

        if delta_method is None:
            delta_method_value = self.delta_method.value
        else:
            delta_method_value = delta_method.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self.num_vol_curves

        # If there is only one time horizon then assume flat vol to this time
        if num_curves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif t_exp <= self.t_exp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif t_exp > self.t_exp[-1]:

            index0 = len(self.t_exp) - 1
            index1 = len(self.t_exp) - 1

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if t_exp <= self.t_exp[i] and t_exp > self.t_exp[i - 1]:
                    index0 = i - 1
                    index1 = i
                    break

        fwd0 = self.fwd[index0]
        fwd1 = self.fwd[index1]

        t0 = self.t_exp[index0]
        t1 = self.t_exp[index1]

        initial_guess = self.k_atm[index0]

        k0 = _solver_for_smile_strike(
            s,
            t_exp,
            self.rd[index0],
            self.rf[index0],
            OptionTypes.EUROPEAN_CALL.value,
            vol_type_value,
            call_delta,
            delta_method_value,
            initial_guess,
            self.parameters[index0],
            self.strikes[index0],
            self.gaps[index0],
        )

        vol0 = vol_function(
            vol_type_value,
            self.parameters[index0],
            self.strikes[index0],
            self.gaps[index0],
            fwd0,
            k0,
            t0,
        )

        if index1 != index0:

            k1 = _solver_for_smile_strike(
                s,
                t_exp,
                self.rd[index1],
                self.rf[index1],
                OptionTypes.EUROPEAN_CALL.value,
                vol_type_value,
                call_delta,
                delta_method_value,
                initial_guess,
                self.parameters[index1],
                self.strikes[index1],
                self.gaps[index1],
            )

            vol1 = vol_function(
                vol_type_value,
                self.parameters[index1],
                self.strikes[index1],
                self.gaps[index1],
                fwd1,
                k1,
                t1,
            )
        else:
            vol1 = vol0

        # In the expiry time dimension, both volatilities are interpolated
        # at the same strikes but different deltas.
        vart0 = vol0 * vol0 * t0
        vart1 = vol1 * vol1 * t1

        if np.abs(t1 - t0) > 1e-6:

            vart = ((t_exp - t0) * vart1 + (t1 - t_exp) * vart0) / (t1 - t0)
            kt = ((t_exp - t0) * k1 + (t1 - t_exp) * k0) / (t1 - t0)

            if vart < 0.0:
                raise FinError(
                    "Failed interpolation due to negative variance."
                )

            volt = np.sqrt(vart / t_exp)

        else:

            volt = vol0
            kt = k0

        return volt, kt

    ###########################################################################

    def _build_vol_surface(
        self, fin_solver_type=FinSolverTypes.NELDER_MEAD, tol=1e-8
    ):
        """Main function to construct the vol surface."""

        s = self.spot_fx_rate
        num_vol_curves = self.num_vol_curves

        if self.vol_func_type == VolFuncTypes.CLARK:
            num_parameters = 3
        elif self.vol_func_type == VolFuncTypes.SABR:
            num_parameters = 4
        elif self.vol_func_type == VolFuncTypes.SABR_BETA_ONE:
            num_parameters = 3
        elif self.vol_func_type == VolFuncTypes.SABR_BETA_HALF:
            num_parameters = 3
        elif self.vol_func_type == VolFuncTypes.BBG:
            num_parameters = 3
        elif self.vol_func_type == VolFuncTypes.CLARK5:
            num_parameters = 5
        else:
            print(self.vol_func_type)
            raise FinError("Unknown Model Type")

        self.parameters = np.zeros([num_vol_curves, num_parameters])

        num_strikes = 5
        self.strikes = np.zeros([num_vol_curves, num_strikes])
        self.gaps = np.zeros([num_vol_curves, num_strikes])

        self.t_exp = np.zeros(num_vol_curves)

        self.fwd = np.zeros(num_vol_curves)
        self.rd = np.zeros(num_vol_curves)
        self.rf = np.zeros(num_vol_curves)
        self.k_atm = np.zeros(num_vol_curves)
        self.delta_atm = np.zeros(num_vol_curves)

        self.k_25d_c = np.zeros(num_vol_curves)
        self.k_25d_p = np.zeros(num_vol_curves)
        self.k_25d_c_ms = np.zeros(num_vol_curves)
        self.k_25d_p_ms = np.zeros(num_vol_curves)
        self.v_25d_ms = np.zeros(num_vol_curves)

        self.k_10d_c = np.zeros(num_vol_curves)
        self.k_10d_p = np.zeros(num_vol_curves)
        self.k_10d_c_ms = np.zeros(num_vol_curves)
        self.k_10d_p_ms = np.zeros(num_vol_curves)
        self.v_10d_ms = np.zeros(num_vol_curves)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        spot_dt = self.value_dt

        for i in range(0, num_vol_curves):

            expiry_dt = self.expiry_dts[i]
            t_exp = (expiry_dt - spot_dt) / g_days_in_year

            dom_df = self.domestic_curve.df(expiry_dt)
            for_df = self.foreign_curve.df(expiry_dt)
            f = s * for_df / dom_df

            self.t_exp[i] = t_exp
            self.rd[i] = -np.log(dom_df) / t_exp
            self.rf[i] = -np.log(for_df) / t_exp
            self.fwd[i] = f

            atm_vol = self.atm_vols[i]

            # This follows exposition in Clarke Page 52
            if self.atm_method == FinFXATMMethod.SPOT:
                self.k_atm[i] = s
            elif self.atm_method == FinFXATMMethod.FWD:
                self.k_atm[i] = f
            elif self.atm_method == FinFXATMMethod.FWD_DELTA_NEUTRAL:
                self.k_atm[i] = f * np.exp(atm_vol * atm_vol * t_exp / 2.0)
            elif self.atm_method == FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ:
                self.k_atm[i] = f * np.exp(-atm_vol * atm_vol * t_exp / 2.0)
            else:
                raise FinError("Unknown Delta Type")

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        ginit = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

        x_inits = []
        ginits = []

        for i in range(0, num_vol_curves):

            atm_vol = self.atm_vols[i]

            if self.use_ms_25d_vol > 0:
                ms25 = self.ms_25_delta_vols[i]
            else:
                ms25 = 0.0

            if self.use_rr_25d_vol > 0:
                rr25 = self.rr_25_delta_vols[i]
            else:
                rr25 = 0.0

            if self.use_ms_10d_vol > 0:
                ms10 = self.ms_10_delta_vols[i]
            else:
                ms10 = 0.0

            if self.use_rr_10d_vol > 0:
                rr10 = self.rr_10_delta_vols[i]
            else:
                rr10 = 0.0

            # https://quantpie.co.uk/fx/fx_rr_str.php
            s25 = atm_vol + ms25 + rr25 / 2.0  # 25d Call
            s50 = atm_vol  # ATM
            s75 = atm_vol + ms25 - rr25 / 2.0  # 25d Put (75D Call)

            s10 = atm_vol + ms10 + rr10 / 2.0  # 10d Call
            s50 = atm_vol  # ATM
            s90 = atm_vol + ms10 - rr10 / 2.0  # 10d Put (90D Call)

            if self.vol_func_type == VolFuncTypes.CLARK:

                # Our preference is to fit to the 10d wings first
                if self.use_ms_10d_vol is False:
                    # Fit to 25d
                    c0 = np.log(atm_vol)
                    c1 = 2.0 * np.log(s75 / s25)
                    c2 = 8.0 * np.log(s25 * s75 / atm_vol / atm_vol)
                    x_init = [c0, c1, c2]
                else:
                    # Fit to 10d
                    c0 = np.log(atm_vol)
                    c1 = np.log(s90 / s10) / 0.80
                    c2 = np.log(s10 * s90 / atm_vol / atm_vol) / 0.32
                    x_init = [c0, c1, c2]

            elif self.vol_func_type == VolFuncTypes.SABR:
                # SABR parameters are alpha, nu, rho
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0
                rho = -0.112
                nu = 0.817

                x_init = [alpha, beta, rho, nu]

            elif self.vol_func_type == VolFuncTypes.SABR_BETA_ONE:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 1.0  # FIXED
                rho = -0.112
                nu = 0.817

                x_init = [alpha, rho, nu]

            elif self.vol_func_type == VolFuncTypes.SABR_BETA_HALF:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                beta = 0.50  # FIXED
                rho = -0.112
                nu = 0.817

                x_init = [alpha, rho, nu]

            elif self.vol_func_type == VolFuncTypes.BBG:

                # Our preference is to fit to the 10d wings first
                if self.use_ms_10d_vol is False:
                    # BBG Params if we fit to 25d
                    a = 8.0 * s75 - 16.0 * s50 + 8.0 * s25
                    b = -6.0 * s75 + 16.0 * s50 - 10.0 * s25
                    c = s75 - 3.0 * s50 + 3.0 * s25
                else:
                    # BBG Params if we fit to 10d
                    a = (25.0 * s90 - 50.0 * s50 + 25.0 * s10) / 8.0
                    b = (-15.0 * s90 + 50.0 * s50 - 35.0 * s10) / 8.0
                    c = (5.0 * s90 - 18.0 * s50 + 45.0 * s10) / 32.0

                x_init = [a, b, c]

            elif self.vol_func_type == VolFuncTypes.CLARK5:

                # Our preference is to fit to the 10d wings first
                if self.use_ms_10d_vol is False:
                    # Fit to 25d
                    c0 = np.log(atm_vol)
                    c1 = 2.0 * np.log(s75 / s25)
                    c2 = 8.0 * np.log(s25 * s75 / atm_vol / atm_vol)
                    x_init = [c0, c1, c2, 0.0, 0.0]
                else:
                    # Fit to 10d
                    c0 = np.log(atm_vol)
                    c1 = np.log(s90 / s10) / 0.80
                    c2 = np.log(s10 * s90 / atm_vol / atm_vol) / 0.32
                    x_init = [c0, c1, c2, 0.0, 0.0]

            else:
                raise FinError("Unknown Model Type")

            x_inits.append(x_init)
            ginits.append(ginit)

        delta_method_value = self.delta_method.value
        vol_type_value = self.vol_func_type.value

        for i in range(0, num_vol_curves):

            t = self.t_exp[i]
            r_d = self.rd[i]
            r_f = self.rf[i]
            k_atm = self.k_atm[i]
            atm_vol = self.atm_vols[i]

            # If the data has not been provided, pass a dummy value
            # as I don't want more arguments and Numpy needs floats
            if self.use_ms_25d_vol:
                ms_25d_vol = self.ms_25_delta_vols[i]
                rr_25d_vol = self.rr_25_delta_vols[i]
            else:
                ms_25d_vol = -999.0
                rr_25d_vol = -999.0

            if self.use_ms_10d_vol:
                ms_10d_vol = self.ms_10_delta_vols[i]
                rr_10d_vol = self.rr_10_delta_vols[i]
            else:
                ms_10d_vol = -999.0
                rr_10d_vol = -999.0

            res = _solve_to_horizon(
                s,
                t,
                r_d,
                r_f,
                k_atm,
                atm_vol,
                ms_25d_vol,
                rr_25d_vol,
                ms_10d_vol,
                rr_10d_vol,
                delta_method_value,
                vol_type_value,
                self.alpha,
                x_inits[i],
                ginits[i],
                fin_solver_type,
                tol,
            )

            (
                self.parameters[i, :],
                self.strikes[i, :],
                self.gaps[i:],
                self.k_25d_c_ms[i],
                self.k_25d_p_ms[i],
                self.k_25d_c[i],
                self.k_25d_p[i],
                self.k_10d_c_ms[i],
                self.k_10d_p_ms[i],
                self.k_10d_c[i],
                self.k_10d_p[i],
            ) = res

    ###########################################################################

    def check_calibration(self, verbose: bool, tol: float = 1e-6):
        """Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals."""

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self.value_dt)
            print("SPOT FX RATE:", self.spot_fx_rate)
            print("ALPHA WEIGHT:", self.alpha)
            print("ATM METHOD:", self.atm_method)
            print("DELTA METHOD:", self.delta_method)
            print("==========================================================")

        k_dummy = 999

        for i in range(0, self.num_vol_curves):

            expiry_dt = self.expiry_dts[i]

            if verbose:
                print("TENOR:", self.tenors[i])
                print("EXPIRY DATE:", expiry_dt)
                print("IN ATM VOL: %9.6f %%" % (100.0 * self.atm_vols[i]))

                if self.use_ms_25d_vol:
                    print(
                        "IN MKT STRANGLE 25d VOL: %9.6f %%"
                        % (100.0 * self.ms_25_delta_vols[i])
                    )
                    print(
                        "IN RSK REVERSAL 25d VOL: %9.6f %%"
                        % (100.0 * self.rr_25_delta_vols[i])
                    )

                if self.use_ms_10d_vol:
                    print(
                        "IN MKT STRANGLE 10d VOL: %9.6f %%"
                        % (100.0 * self.ms_10_delta_vols[i])
                    )
                    print(
                        "IN RSK REVERSAL 10d VOL: %9.6f %%"
                        % (100.0 * self.rr_10_delta_vols[i])
                    )

            call = FXVanillaOption(
                expiry_dt,
                k_dummy,
                self.currency_pair,
                OptionTypes.EUROPEAN_CALL,
                1.0,
                self.notional_currency,
            )

            put = FXVanillaOption(
                expiry_dt,
                k_dummy,
                self.currency_pair,
                OptionTypes.EUROPEAN_PUT,
                1.0,
                self.notional_currency,
            )

            ###################################################################
            # AT THE MONEY
            ###################################################################

            if verbose:
                print("======================================================")
                print("T_(YEARS): ", self.t_exp[i])
                print("CNT_cPD_RD:%9.6f %%" % (self.rd[i] * 100))
                print("CNT_cPD_RF:%9.6f %%" % (self.rf[i] * 100))
                print("FWD_RATE:  %9.6f" % (self.fwd[i]))

            sigma_atm_out = vol_function(
                self.vol_func_type.value,
                self.parameters[i],
                self.strikes[i],
                self.gaps[i],
                self.fwd[i],
                self.k_atm[i],
                self.t_exp[i],
            )

            if verbose:
                print("======================================================")
                print("VOL FUNCTION", self.vol_func_type)
                print("VOL_pARAMETERS:", self.parameters[i])
                print("======================================================")
                print("OUT_k_atm:  %9.6f" % (self.k_atm[i]))
                print("OUT_ATM_VOL: %9.6f %%" % (100.0 * sigma_atm_out))

            diff = sigma_atm_out - self.atm_vols[i]

            if np.abs(diff) > tol:
                print(
                    "FAILED FIT TO ATM VOL IN: %9.6f  OUT: %9.6f  DIFF: %9.6f"
                    % (
                        self.atm_vols[i] * 100.0,
                        sigma_atm_out * 100.0,
                        diff * 100.0,
                    )
                )

            call.strike_fx_rate = self.k_atm[i]
            put.strike_fx_rate = self.k_atm[i]

            model = BlackScholes(sigma_atm_out)

            delta_call = call.delta(
                self.value_dt,
                self.spot_fx_rate,
                self.domestic_curve,
                self.foreign_curve,
                model,
            )[self.delta_method_string]

            delta_put = put.delta(
                self.value_dt,
                self.spot_fx_rate,
                self.domestic_curve,
                self.foreign_curve,
                model,
            )[self.delta_method_string]

            if verbose:
                print(
                    "CALL_DELTA: % 9.6f  PUT_DELTA: % 9.6f  NET_DELTA: % 9.6f"
                    % (delta_call, delta_put, delta_call + delta_put)
                )

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.25/-0.25
            ###################################################################

            if self.use_ms_25d_vol is True:

                ms_vol = self.atm_vols[i] + self.ms_25_delta_vols[i]

                if verbose:

                    print("==================================================")
                    print(
                        "MKT STRANGLE 25d VOL IN: %9.6f %%"
                        % (100.0 * self.ms_25_delta_vols[i])
                    )

                call.strike_fx_rate = self.k_25d_c_ms[i]
                put.strike_fx_rate = self.k_25d_p_ms[i]

                model = BlackScholes(ms_vol)

                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                if verbose:
                    print(
                        "k_25d_c_ms: %9.6f  ATM + ms_vol: %9.6f %%   DELTA: %9.6f"
                        % (self.k_25d_c_ms[i], 100.0 * ms_vol, delta_call)
                    )

                    print(
                        "k_25d_p_ms: %9.6f  ATM + ms_vol: %9.6f %%   DELTA: %9.6f"
                        % (self.k_25d_p_ms[i], 100.0 * ms_vol, delta_put)
                    )

                call_value = call.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                put_value = put.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                mkt_strangle_value = call_value + put_value

                if verbose:
                    print(
                        "CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_VALUE: % 9.6f"
                        % (call_value, put_value, mkt_strangle_value)
                    )

                ###############################################################
                # NOW WE ASSIGN A DIFFERENT VOLATILITY TO THE MS STRIKES
                # THE DELTAS WILL NO LONGER EQUAL 0.25, -0.25
                ###############################################################

                # CALL
                sigma_k_25d_c_ms = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_25d_c_ms[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_25d_c_ms)
                call_value = call.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                # THIS IS NOT GOING TO BE 0.25 AS WE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                # PUT
                sigma_k_25d_p_ms = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_25d_p_ms[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_25d_p_ms)
                put_value = put.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                # THIS IS NOT GOING TO BE -0.25 AS WE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                mkt_strangle_value_skew = call_value + put_value

                if verbose:
                    print(
                        "k_25d_c_ms: %9.6f  SURFACE_VOL: %9.6f %% DELTA: %9.6f"
                        % (
                            self.k_25d_c_ms[i],
                            100.0 * sigma_k_25d_c_ms,
                            delta_call,
                        )
                    )

                    print(
                        "k_25d_p_ms: %9.6f  SURFACE_VOL: %9.6f %% DELTA: %9.6f"
                        % (
                            self.k_25d_p_ms[i],
                            100.0 * sigma_k_25d_p_ms,
                            delta_put,
                        )
                    )

                    print(
                        "CALL_VALUE: %9.6f  PUT_VALUE: %9.6f MS_SKEW_VALUE: % 9.6f"
                        % (call_value, put_value, mkt_strangle_value_skew)
                    )

                diff = mkt_strangle_value - mkt_strangle_value_skew
                if np.abs(diff) > tol:
                    print(
                        "FAILED FIT TO 25d MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f"
                        % (mkt_strangle_value, mkt_strangle_value_skew, diff)
                    )

                ###############################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.25, -0.25
                ###############################################################

                call.strike_fx_rate = self.k_25d_c[i]
                put.strike_fx_rate = self.k_25d_p[i]

                sigma_k_25d_c = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_25d_c[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_25d_c)

                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                sigma_k_25d_p = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_25d_p[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_25d_p)

                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                if verbose:
                    print(
                        "k_25d_c: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                        % (self.k_25d_c[i], 100.0 * sigma_k_25d_c, delta_call)
                    )

                    print(
                        "k_25d_p: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                        % (self.k_25d_p[i], 100.0 * sigma_k_25d_p, delta_put)
                    )

                sigma_rr = sigma_k_25d_c - sigma_k_25d_p

                if verbose:
                    print(
                        "=========================================================="
                    )
                    print(
                        "RR = VOL_k_25_c - VOL_k_25_p => RR_IN: %9.6f %% RR_OUT: %9.6f %%"
                        % (100.0 * self.rr_25_delta_vols[i], 100.0 * sigma_rr)
                    )
                    print(
                        "=========================================================="
                    )

                diff = sigma_rr - self.rr_25_delta_vols[i]

                if np.abs(diff) > tol:
                    print(
                        "FAILED FIT TO 25d RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f"
                        % (
                            self.rr_25_delta_vols[i] * 100.0,
                            sigma_rr * 100.0,
                            diff * 100.0,
                        )
                    )

            ###################################################################
            # NOW WE ASSIGN THE SAME VOLATILITY TO THE MS STRIKES
            # THESE STRIKES ARE DETERMINED BY SETTING DELTA TO 0.10/-0.10
            ###################################################################

            if self.use_ms_10d_vol:

                ms_vol = self.atm_vols[i] + self.ms_10_delta_vols[i]

                if verbose:

                    print(
                        "=========================================================="
                    )
                    print(
                        "MKT STRANGLE 10d VOL IN: %9.6f %%"
                        % (100.0 * self.ms_10_delta_vols[i])
                    )

                call.strike_fx_rate = self.k_10d_c_ms[i]
                put.strike_fx_rate = self.k_10d_p_ms[i]

                model = BlackScholes(ms_vol)

                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                if verbose:
                    print(
                        "k_10d_c_ms: %9.6f  ATM + ms_vol: %9.6f %%   DELTA: %9.6f"
                        % (self.k_10d_c_ms[i], 100.0 * ms_vol, delta_call)
                    )

                    print(
                        "k_10d_p_ms: %9.6f  ATM + ms_vol: %9.6f %%   DELTA: %9.6f"
                        % (self.k_10d_p_ms[i], 100.0 * ms_vol, delta_put)
                    )

                call_value = call.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                put_value = put.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                mkt_strangle_value = call_value + put_value

                if verbose:
                    print(
                        "CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_VALUE: % 9.6f"
                        % (call_value, put_value, mkt_strangle_value)
                    )

                ###############################################################
                # NOW WE ASSIGN A DIFFERENT VOLATILITY TO THE MS STRIKES
                # THE DELTAS WILL NO LONGER EQUAL 0.25, -0.25
                ###############################################################

                # CALL
                sigma_k_10d_c_ms = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_10d_c_ms[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_10d_c_ms)
                call_value = call.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                # THIS IS NOT GOING TO BE 0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                # PUT
                sigma_k_10d_p_ms = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_10d_p_ms[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_10d_p_ms)
                put_value = put.value(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )["v"]

                # THIS IS NOT GOING TO BE -0.10 AS WE HAVE USED A DIFFERENT SKEW VOL
                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                mkt_strangle_value_skew = call_value + put_value

                if verbose:
                    print(
                        "k_10d_c_ms: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                        % (
                            self.k_10d_c_ms[i],
                            100.0 * sigma_k_10d_c_ms,
                            delta_call,
                        )
                    )

                    print(
                        "k_10d_p_ms: %9.6f  SURFACE_VOL: %9.6f %%   DELTA: %9.6f"
                        % (
                            self.k_10d_p_ms[i],
                            100.0 * sigma_k_10d_p_ms,
                            delta_put,
                        )
                    )

                    print(
                        "CALL_VALUE: %9.6f  PUT_VALUE: %9.6f  MS_SKEW_VALUE: % 9.6f"
                        % (call_value, put_value, mkt_strangle_value_skew)
                    )

                diff = mkt_strangle_value - mkt_strangle_value_skew
                if np.abs(diff) > tol:
                    print(
                        "FAILED FIT TO 10d MS VAL: %9.6f  OUT: %9.6f  DIFF: % 9.6f"
                        % (mkt_strangle_value, mkt_strangle_value_skew, diff)
                    )

                ###############################################################
                # NOW WE SHIFT STRIKES SO THAT DELTAS NOW EQUAL 0.10, -0.10
                ###############################################################

                call.strike_fx_rate = self.k_10d_c[i]
                put.strike_fx_rate = self.k_10d_p[i]

                sigma_k_10d_c = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_10d_c[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_10d_c)

                # THIS DELTA SHOULD BE +0.25
                delta_call = call.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                sigma_k_10d_p = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i],
                    self.strikes[i],
                    self.gaps[i],
                    self.fwd[i],
                    self.k_10d_p[i],
                    self.t_exp[i],
                )

                model = BlackScholes(sigma_k_10d_p)

                # THIS DELTA SHOULD BE -0.25
                delta_put = put.delta(
                    self.value_dt,
                    self.spot_fx_rate,
                    self.domestic_curve,
                    self.foreign_curve,
                    model,
                )[self.delta_method_string]

                if verbose:
                    print(
                        "k_10d_c: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                        % (self.k_10d_c[i], 100.0 * sigma_k_10d_c, delta_call)
                    )

                    print(
                        "k_10d_p: %9.7f  VOL: %9.6f  DELTA: % 9.6f"
                        % (self.k_10d_p[i], 100.0 * sigma_k_10d_p, delta_put)
                    )

                sigma_rr = sigma_k_10d_c - sigma_k_10d_p

                if verbose:
                    print(
                        "====================================================="
                    )
                    print(
                        "RR = VOL_k_10d_c - VOL_k_10d_p => RR_IN: %9.6f %% RR_OUT: %9.6f %%"
                        % (100.0 * self.rr_10_delta_vols[i], 100.0 * sigma_rr)
                    )
                    print(
                        "====================================================="
                    )

                diff = sigma_rr - self.rr_10_delta_vols[i]

                if np.abs(diff) > tol:
                    print(
                        "FAILED FIT TO 10d RRV IN: % 9.6f  OUT: % 9.6f  DIFF: % 9.6f"
                        % (
                            self.rr_10_delta_vols[i] * 100.0,
                            sigma_rr * 100.0,
                            diff * 100.0,
                        )
                    )

    ###########################################################################

    def implied_dbns(self, low_fx, high_fx, num_intervals):
        """Calculate the pdf for each tenor horizon. Returns a list of
        FinDistribution objects, one for each tenor horizon."""

        dbns = []

        for i_tenor in range(0, len(self.tenors)):

            f = self.fwd[i_tenor]
            t = self.t_exp[i_tenor]

            d_fx = (high_fx - low_fx) / num_intervals

            dom_df = self.domestic_curve.df_t(t)
            for_df = self.foreign_curve.df_t(t)

            r_d = -np.log(dom_df) / t
            r_f = -np.log(for_df) / t

            ks = []
            vols = []

            for ik in range(0, num_intervals):

                k = low_fx + ik * d_fx

                vol = vol_function(
                    self.vol_func_type.value,
                    self.parameters[i_tenor],
                    self.strikes[i_tenor],
                    self.gaps[i_tenor],
                    f,
                    k,
                    t,
                )

                ks.append(k)
                vols.append(vol)

            ks = np.array(ks)
            vols = np.array(vols)

            density = option_implied_dbn(
                self.spot_fx_rate, t, r_d, r_f, ks, vols
            )

            dbn = FinDistribution(ks, density)
            dbns.append(dbn)

        return dbns

    ###########################################################################

    def plot_vol_curves(self):
        """Generates a plot of each of the vol discount implied by the market
        and fitted."""

        plt.figure()

        vol_type_val = self.vol_func_type.value

        for tenor_index in range(0, self.num_vol_curves):

            atm_vol = self.atm_vols[tenor_index] * 100
            ms_vol25 = self.ms_25_delta_vols[tenor_index] * 100
            rr_vol_25 = self.rr_25_delta_vols[tenor_index] * 100
            ms_vol10 = self.ms_10_delta_vols[tenor_index] * 100
            rr_vol_10 = self.rr_10_delta_vols[tenor_index] * 100
            strikes = self.strikes[tenor_index]

            gaps = self.gaps[tenor_index]

            low_k = self.k_10d_p[tenor_index] * 0.90
            high_k = self.k_10d_c_ms[tenor_index] * 1.10

            ks = []
            vols = []
            num_intervals = 30
            k = low_k
            dk = (high_k - low_k) / num_intervals
            params = self.parameters[tenor_index]
            t = self.t_exp[tenor_index]
            f = self.fwd[tenor_index]

            for _ in range(0, num_intervals):

                sigma = (
                    vol_function(vol_type_val, params, strikes, gaps, f, k, t)
                    * 100.0
                )
                ks.append(k)
                vols.append(sigma)
                k = k + dk

            label_str = self.tenors[tenor_index]
            label_str += " ATM: " + str(atm_vol)[0:6]
            label_str += " MS25: " + str(ms_vol25)[0:6]
            label_str += " RR25: " + str(rr_vol_25)[0:6]
            label_str += " MS10: " + str(ms_vol10)[0:6]
            label_str += " RR10: " + str(rr_vol_10)[0:6]

            plt.plot(ks, vols, label=label_str)
            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = (
                "JNT FIT:" + self.currency_pair + " " + str(self.vol_func_type)
            )

            key_strikes = []
            key_strikes.append(self.k_atm[tenor_index])

            key_vols = []
            for K in key_strikes:

                sigma = (
                    vol_function(vol_type_val, params, strikes, gaps, f, K, t)
                    * 100.0
                )

                key_vols.append(sigma)

            plt.plot(key_strikes, key_vols, "ko", markersize=4)

            key_strikes = []
            key_strikes.append(self.k_25d_p[tenor_index])
            key_strikes.append(self.k_25d_p_ms[tenor_index])
            key_strikes.append(self.k_25d_c[tenor_index])
            key_strikes.append(self.k_25d_c_ms[tenor_index])

            key_vols = []
            for K in key_strikes:

                sigma = (
                    vol_function(vol_type_val, params, strikes, gaps, f, K, t)
                    * 100.0
                )

                key_vols.append(sigma)

            plt.plot(key_strikes, key_vols, "bo", markersize=4)

            key_strikes = []
            key_strikes.append(self.k_10d_p[tenor_index])
            key_strikes.append(self.k_10d_p_ms[tenor_index])
            key_strikes.append(self.k_10d_c[tenor_index])
            key_strikes.append(self.k_10d_c_ms[tenor_index])

            key_vols = []
            for K in key_strikes:
                sigma = (
                    vol_function(vol_type_val, params, strikes, gaps, f, K, t)
                    * 100.0
                )
                key_vols.append(sigma)

            plt.plot(key_strikes, key_vols, "ro", markersize=4)

        plt.title(title)
        plt.legend(loc="lower left", bbox_to_anchor=(1, 0))

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUE DATE", self.value_dt)
        s += label_to_string("FX RATE", self.spot_fx_rate)
        s += label_to_string("CCY PAIR", self.currency_pair)
        s += label_to_string("NOTIONAL CCY", self.notional_currency)
        s += label_to_string("NUM TENORS", self.num_vol_curves)
        s += label_to_string("ATM METHOD", self.atm_method)
        s += label_to_string("DELTA METHOD", self.delta_method)
        s += label_to_string("ALPHA WEIGHT", self.alpha)
        s += label_to_string("VOL FUNCTION", self.vol_func_type)

        for i in range(0, self.num_vol_curves):

            s += "\n"

            s += label_to_string("TENOR", self.tenors[i])
            s += label_to_string("EXPIRY DATE", self.expiry_dts[i])
            s += label_to_string("TIME (YRS)", self.t_exp[i])
            s += label_to_string("FWD FX", self.fwd[i])

            s += label_to_string("ATM VOLS", self.atm_vols[i] * 100.0)
            s += label_to_string("MS VOLS", self.ms_25_delta_vols[i] * 100.0)
            s += label_to_string("RR VOLS", self.rr_25_delta_vols[i] * 100.0)

            s += label_to_string("ATM Strike", self.k_atm[i])
            s += label_to_string("ATM Delta", self.delta_atm[i])

            s += label_to_string("k_atm", self.k_atm[i])

            s += label_to_string("MS 25d Call Strike", self.k_25d_c_ms[i])
            s += label_to_string("MS 25d Put Strike", self.k_25d_p_ms[i])
            s += label_to_string("SKEW 25d CALL STRIKE", self.k_25d_c[i])
            s += label_to_string("SKEW 25d PUT STRIKE", self.k_25d_p[i])
            s += label_to_string("PARAMS", self.parameters[i])

            s += label_to_string("MS 10d Call Strike", self.k_10d_c_ms[i])
            s += label_to_string("MS 10d Put Strike", self.k_10d_p_ms[i])
            s += label_to_string("SKEW 10d CALL STRIKE", self.k_10d_c[i])
            s += label_to_string("SKEW 10d PUT STRIKE", self.k_10d_p[i])

        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
