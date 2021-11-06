##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane, Saeed Amen
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

###############################################################################
# TODO: Speed up search for strike by providing derivative function to go with
#       delta fit.
###############################################################################


@njit(fastmath=True, cache=True)
def g(K, *args):
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
# Do not cache this function


@njit(fastmath=True)  # , cache=True)
def obj_fast(params, *args):
    """ Return a function that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by cvec
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
    delta_method_value = args[9]
    targetRRVol = args[10]
    vol_type_value = args[11]

    f = s * np.exp((rd-rf)*t)
    # We first need to solve for the strikes at the 25 delta points using the
    # new volatility curve

    # Match the at-the-money option volatility
    atm_curve_vol = vol_function(vol_type_value, params, f, K_ATM, t)
    term1 = (atm_vol - atm_curve_vol)**2

    ###########################################################################
    # Match the market strangle value but this has to be at the MS strikes
    ###########################################################################

    sigma_K_25D_C_MS = vol_function(vol_type_value, params, f, K_25D_C_MS, t)

    V_25D_C_MS = bs_value(s, t, K_25D_C_MS, rd, rf, sigma_K_25D_C_MS,
                          OptionTypes.EUROPEAN_CALL.value)

    sigma_K_25D_P_MS = vol_function(vol_type_value, params, f, K_25D_P_MS, t)

    V_25D_P_MS = bs_value(s, t, K_25D_P_MS, rd, rf, sigma_K_25D_P_MS,
                          OptionTypes.EUROPEAN_PUT.value)

    V_25D_MS = V_25D_C_MS + V_25D_P_MS
    term2 = (V_25D_MS - V_25D_MS_target)**2

    ###########################################################################
    # Match the risk reversal volatility
    ###########################################################################

    K_25D_C = solver_for_smile_strike_fast(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_CALL.value,
                                           vol_type_value, +0.2500,
                                           delta_method_value, K_25D_C_MS,
                                           params)

    sigma_K_25D_C = vol_function(vol_type_value, params, f, K_25D_C, t)

    K_25D_P = solver_for_smile_strike_fast(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_PUT.value,
                                           vol_type_value, -0.2500,
                                           delta_method_value, K_25D_P_MS,
                                           params)

    sigma_K_25D_P = vol_function(vol_type_value, params, f, K_25D_P, t)

    sigma_25D_RR = (sigma_K_25D_C - sigma_K_25D_P)
    term3 = (sigma_25D_RR - targetRRVol)**2

    # sum up the errors
    err = term1 + term2 + term3

    return err

###############################################################################
# This function cannot be jitted until the scipy minimisation has been replaced
# with a jittable function


def solve_to_horizon_fast(s, t,
                          rd, rf,
                          K_ATM, atm_vol,
                          ms25DVol, rr25DVol,
                          delta_method_value, vol_type_value,
                          xopt):

    c0 = xopt

    # Determine the price of a market strangle from market strangle
    # Need to price a call and put that agree with market strangle

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

    # Determine parameters of vol surface using minimisation
    tol = 1e-8

    fargs = (s, t, rd, rf,
             K_ATM, atm_vol,
             K_25D_C_MS, K_25D_P_MS,
             V_25D_MS,
             delta_method_value, rr25DVol, vol_type_value)

    opt = minimize(obj_fast, c0, fargs, method="CG", tol=tol)
    xopt = opt.x

    params = np.array(xopt)

    K_25D_C = solver_for_smile_strike_fast(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_CALL.value,
                                           vol_type_value, +0.2500,
                                           delta_method_value, K_25D_C_MS,
                                           params)

    K_25D_P = solver_for_smile_strike_fast(s, t, rd, rf,
                                           OptionTypes.EUROPEAN_PUT.value,
                                           vol_type_value, -0.2500,
                                           delta_method_value, K_25D_P_MS,
                                           params)

    ret = (params, K_25D_C_MS, K_25D_P_MS, K_25D_C, K_25D_P)
    return ret

###############################################################################


@njit(float64(int64, float64[:], float64, float64, float64),
      cache=True, fastmath=True)
def vol_function(vol_function_type_value, params, f, k, t):
    """ Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book. """

    if vol_function_type_value == VolFunctionTypes.CLARK.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR.value:
        vol = vol_function_sabr(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_ONE.value:
        vol = vol_function_sabr_beta_one(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.SABR_BETA_HALF.value:
        vol = vol_function_sabr_beta_half(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.BBG.value:
        vol = vol_function_bloomberg(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFunctionTypes.CLARK5.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")

###############################################################################


@njit(cache=True, fastmath=True)
def delta_fit(K, *args):
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

    f = s*np.exp((rd-rf)*t)
    v = vol_function(vol_type_value, params, f, K, t)
    delta_out = fast_delta(
        s, t, K, rd, rf, v, deltaTypeValue, option_type_value)
    inverseDeltaOut = norminvcdf(np.abs(delta_out))
    invObjFn = inverseDeltaTarget - inverseDeltaOut

    return invObjFn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.


@njit(float64(float64, float64, float64, float64, int64, int64, float64,
              int64, float64, float64[:]), fastmath=True)
def solver_for_smile_strike_fast(s, t, rd, rf,
                                 option_type_value,
                                 volatilityTypeValue,
                                 delta_target,
                                 delta_method_value,
                                 initialGuess,
                                 parameters):
    """ Solve for the strike that sets the delta of the option equal to the
    target value of delta allowing the volatility to be a function of the
    strike. """

    inverseDeltaTarget = norminvcdf(np.abs(delta_target))

    argtuple = (volatilityTypeValue, s, t, rd, rf, option_type_value,
                delta_method_value, inverseDeltaTarget, parameters)

    K = newton_secant(delta_fit, x0=initialGuess, args=argtuple,
                      tol=1e-8, maxiter=50)

    return K

###############################################################################
# Unable to cache function


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
        K = F0T * np.exp(-vsqrtt * (phi*norminvdelta - vsqrtt/2.0))
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
        arg = delta_target*phi   # CHECK THIS!!!!!!!!
        norminvdelta = norminvcdf(arg)
        K = F0T * np.exp(-vsqrtt * (phi*norminvdelta - vsqrtt/2.0))
        return K

    elif delta_method_value == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

        argtuple = (spot_fx_rate, tdel, rd, rf, volatility,
                    delta_method_value, option_type_value, delta_target)

        K = newton_secant(g, x0=spot_fx_rate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    elif delta_method_value == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

        argtuple = (spot_fx_rate, tdel, rd, rf, volatility,
                    delta_method_value, option_type_value, delta_target)

        K = newton_secant(g, x0=spot_fx_rate, args=argtuple,
                          tol=1e-7, maxiter=50)

        return K

    else:

        raise FinError("Unknown FinFXDeltaMethod")

###############################################################################


class FXVolSurface():
    """ Class to perform a calibration of a chosen parametrised surface to the
    prices of FX options at different strikes and expiry tenors. The 
    calibration inputs are the ATM and 25 Delta volatilities given in terms of
    the market strangle amd risk reversals. There is a choice of volatility
    function ranging from polynomial in delta to a limited version of SABR. """

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
                 atmMethod: FinFXATMMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL,
                 deltaMethod: FinFXDeltaMethod = FinFXDeltaMethod.SPOT_DELTA,
                 volatility_function_type: VolFunctionTypes = VolFunctionTypes.CLARK):
        """ Create the FinFXVolSurface object by passing in market vol data
        for ATM and 25 Delta Market Strangles and Risk Reversals. """

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

        if len(atm_vols) != self._num_vol_curves:
            raise FinError("Number ATM vols must equal number of tenors")

        if len(atm_vols) != self._num_vol_curves:
            raise FinError("Number ATM vols must equal number of tenors")

        if len(mktStrangle25DeltaVols) != self._num_vol_curves:
            raise FinError("Number MS25D vols must equal number of tenors")

        if len(riskReversal25DeltaVols) != self._num_vol_curves:
            raise FinError("Number RR25D vols must equal number of tenors")

        self._tenors = tenors
        self._atm_vols = np.array(atm_vols)/100.0
        self._mktStrangle25DeltaVols = np.array(mktStrangle25DeltaVols)/100.0
        self._riskReversal25DeltaVols = np.array(riskReversal25DeltaVols)/100.0

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

        self.build_vol_surface()

###############################################################################

    def volatility(self, K, expiry_date):
        """ Interpolate the Black-Scholes volatility from the volatility
        surface given the option strike and expiry date. Linear interpolation
        is done in variance x time. """

        vol_type_value = self._volatility_function_type.value

        index0 = 0
        index1 = 0

        t = (expiry_date - self._valuation_date) / gDaysInYear

        num_curves = self._num_vol_curves

        if num_curves == 1:

            # The volatility term structure is flat if there is only one expiry
            fwd = self._F0T[0]
            texp = self._texp[0]
            vol = vol_function(vol_type_value, self._parameters[0],
                               fwd, K, texp)
            return vol

        # If the time is below first time then assume a flat vol
        if t <= self._texp[0]:

            fwd = self._F0T[0]
            texp = self._texp[0]
            vol = vol_function(vol_type_value, self._parameters[0],
                               fwd, K, texp)
            return vol

        # If the time is beyond the last time then extrapolate with a flat vol
        if t > self._texp[-1]:

            fwd = self._F0T[-1]
            texp = self._texp[-1]
            vol = vol_function(vol_type_value, self._parameters[-1],
                               fwd, K, texp)
            return vol

        for i in range(1, num_curves):

            if t <= self._texp[i] and t > self._texp[i-1]:
                index0 = i-1
                index1 = i
                break

        fwd0 = self._F0T[index0]
        t0 = self._texp[index0]
        vol0 = vol_function(vol_type_value, self._parameters[index0],
                            fwd0, K, t0)

        fwd1 = self._F0T[index1]
        t1 = self._texp[index1]
        vol1 = vol_function(vol_type_value, self._parameters[index1],
                            fwd1, K, t1)

        vart0 = vol0*vol0*t0
        vart1 = vol1*vol1*t1
        vart = ((t-t0) * vart1 + (t1-t) * vart0) / (t1 - t0)

        if vart < 0.0:
            raise FinError("Negative variance.")

        volt = np.sqrt(vart/t)
        return volt

###############################################################################

    def build_vol_surface(self):

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
        self._texp = np.zeros(num_vol_curves)

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

        x_inits = []
        for i in range(0, num_vol_curves):

            atm_vol = self._atm_vols[i]
            ms25 = self._mktStrangle25DeltaVols[i]
            rr25 = self._riskReversal25DeltaVols[i]
            s25 = atm_vol + ms25 + rr25/2.0
            s50 = atm_vol
            s75 = atm_vol + ms25 - rr25/2.0

            if self._volatility_function_type == VolFunctionTypes.CLARK:

                # Fit to 25D
                c0 = np.log(atm_vol)
                c1 = 2.0 * np.log(s75/s25)
                c2 = 8.0 * np.log(s25*s75/atm_vol/atm_vol)
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
                rho = -0.112
                nu = 0.817
                x_init = [alpha, nu, rho]

            elif self._volatility_function_type == VolFunctionTypes.SABR_BETA_HALF:
                # SABR parameters are alpha, nu, rho
                alpha = 0.174
                rho = -0.112
                nu = 0.817
                x_init = [alpha, rho, nu]

            elif self._volatility_function_type == VolFunctionTypes.BBG:

                # BBG Params if we fit to 25D
                a = 8.0*s75-16.0*s50+8.0*s25
                b = -6.0*s75+16.0*s50-10.0*s25
                c = s75-3.0*s50+3.0*s25

                x_init = [a, b, c]

            elif self._volatility_function_type == VolFunctionTypes.CLARK5:

                # Fit to 25D
                c0 = np.log(atm_vol)
                c1 = 2.0 * np.log(s75/s25)
                c2 = 8.0 * np.log(s25*s75/atm_vol/atm_vol)
                x_init = [c0, c1, c2, 0.0, 0.0]

            else:
                raise FinError("Unknown Model Type")

            x_inits.append(x_init)

        delta_method_value = self._deltaMethod.value
        vol_type_value = self._volatility_function_type.value

        for i in range(0, num_vol_curves):

            t = self._texp[i]
            rd = self._rd[i]
            rf = self._rf[i]
            K_ATM = self._K_ATM[i]
            atm_vol = self._atm_vols[i]
            ms25DVol = self._mktStrangle25DeltaVols[i]
            rr25DVol = self._riskReversal25DeltaVols[i]

#            print(t, rd, rf, K_ATM, atm_vol, ms25DVol, rr25DVol)

            res = solve_to_horizon_fast(s, t, rd, rf, K_ATM,
                                        atm_vol, ms25DVol, rr25DVol,
                                        delta_method_value, vol_type_value,
                                        x_inits[i])

            (self._parameters[i, :],
             self._K_25D_C_MS[i], self._K_25D_P_MS[i],
             self._K_25D_C[i], self._K_25D_P[i]) = res

###############################################################################

    def solver_for_smile_strike(self,
                                option_type_value,
                                delta_target,
                                tenorIndex,
                                initialValue):
        """ Solve for the strike that sets the delta of the option equal to the
        target value of delta allowing the volatility to be a function of the
        strike. """

        s0 = self._spot_fx_rate
        tdel = self._texp[tenorIndex]
        rd = self._rd[tenorIndex]
        rf = self._rf[tenorIndex]

        inverseDeltaTarget = norminvcdf(np.abs(delta_target))
        argtuple = (self, s0, tdel, rd, rf, option_type_value,
                    inverseDeltaTarget, tenorIndex)

        vol_type_value = self._volatility_function_type.value

        argtuple = (vol_type_value, s0, tdel, rd, rf, option_type_value,
                    self._deltaMethod.value,
                    inverseDeltaTarget, self._parameters[tenorIndex])

        K = newton_secant(delta_fit, x0=initialValue, args=argtuple,
                          tol=1e-5, maxiter=50)

        return K

###############################################################################

    def check_calibration(self, verbose: bool, tol: float = 1e-6):

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self._valuation_date)
            print("SPOT FX RATE:", self._spot_fx_rate)
            print("ATM METHOD:", self._atmMethod)
            print("DELTA METHOD:", self._deltaMethod)
            print("==========================================================")

        K_dummy = 999

        for i in range(0, self._num_vol_curves):

            expiry_date = self._expiry_dates[i]

            if verbose:
                print("TENOR:", self._tenors[i])
                print("EXPIRY DATE:", expiry_date)
                print("IN ATM VOL: %9.6f %%" % (100.0*self._atm_vols[i]))
                print("IN MKT STRANGLE 25D VOL: %9.6f %%" %
                      (100.0*self._mktStrangle25DeltaVols[i]))
                print("IN RSK REVERSAL 25D VOL: %9.6f %%" %
                      (100.0*self._riskReversal25DeltaVols[i]))

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

            msVol = self._atm_vols[i] + self._mktStrangle25DeltaVols[i]

            if verbose:

                print("==========================================================")
                print("MKT STRANGLE VOL IN: %9.6f %%"
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

###############################################################################

    def implied_dbns(self, lowFX, highFX, numIntervals):
        """ Calculate the pdf for each tenor horizon. Returns a list of 
        FinDistribution objects, one for each tenor horizon. """

        dbns = []

        for iTenor in range(0, len(self._tenors)):

            f = self._F0T[iTenor]
            texp = self._texp[iTenor]

            dFX = (highFX - lowFX) / numIntervals

            domDF = self._dom_discount_curve._df(texp)
            forDF = self._for_discount_curve._df(texp)

            rd = -np.log(domDF) / texp
            rf = -np.log(forDF) / texp

            Ks = []
            vols = []

            for iK in range(0, numIntervals):

                k = lowFX + iK*dFX

                vol = vol_function(self._volatility_function_type.value,
                                   self._parameters[iTenor],
                                   f, k, texp)

                Ks.append(k)
                vols.append(vol)

            Ks = np.array(Ks)
            vols = np.array(vols)

            density = option_implied_dbn(self._spot_fx_rate, texp,
                                         rd, rf, Ks, vols)

            dbn = FinDistribution(Ks, density)
            dbns.append(dbn)

        return dbns

###############################################################################

    def plot_vol_curves(self):

        plt.figure()

        volTypeVal = self._volatility_function_type.value

        for tenorIndex in range(0, self._num_vol_curves):

            atm_vol = self._atm_vols[tenorIndex]*100
            msVol = self._mktStrangle25DeltaVols[tenorIndex]*100
            rrVol = self._riskReversal25DeltaVols[tenorIndex]*100

            lowK = self._K_25D_P[tenorIndex] * 0.75
            highK = self._K_25D_C[tenorIndex] * 1.25

            strikes = []
            vols = []
            numIntervals = 30
            K = lowK
            dK = (highK - lowK)/numIntervals
            params = self._parameters[tenorIndex]
            t = self._texp[tenorIndex]
            f = self._F0T[tenorIndex]

            for _ in range(0, numIntervals):
                sigma = vol_function(volTypeVal, params, f, K, t) * 100.0
                strikes.append(K)
                vols.append(sigma)
                K = K + dK

            labelStr = self._tenors[tenorIndex]
            labelStr += " ATM: " + str(atm_vol)[0:6]
            labelStr += " MS: " + str(msVol)[0:6]
            labelStr += " RR: " + str(rrVol)[0:6]

            plt.plot(strikes, vols, label=labelStr)
            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = "25D FIT:" + self._currency_pair + \
                " " + str(self._volatility_function_type)

            keyStrikes = []
            keyStrikes.append(self._K_ATM[tenorIndex])

            keyVols = []
            for K in keyStrikes:
                sigma = vol_function(volTypeVal, params, f, K, t) * 100.0
                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'ko', markersize=4)

            keyStrikes = []
            keyStrikes.append(self._K_25D_P[tenorIndex])
            keyStrikes.append(self._K_25D_P_MS[tenorIndex])
            keyStrikes.append(self._K_25D_C[tenorIndex])
            keyStrikes.append(self._K_25D_C_MS[tenorIndex])

            keyVols = []

            for K in keyStrikes:
                sigma = vol_function(volTypeVal, params, f, K, t) * 100.0
                keyVols.append(sigma)

            plt.plot(keyStrikes, keyVols, 'bo', markersize=4)

        plt.title(title)
#        plt.legend(loc="lower left", bbox_to_anchor=(1,0))

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
        s += label_to_string("VOL FUNCTION", self._volatility_function_type)

        for i in range(0, self._num_vol_curves):

            s += "\n"

            s += label_to_string("TENOR", self._tenors[i])
            s += label_to_string("EXPIRY DATE", self._expiry_dates[i])
            s += label_to_string("TIME (YRS)", self._texp[i])
            s += label_to_string("FWD FX", self._F0T[i])

            s += label_to_string("ATM VOLS", self._atm_vols[i]*100.0)
            s += label_to_string("MS VOLS",
                                 self._mktStrangle25DeltaVols[i]*100.0)
            s += label_to_string("RR VOLS",
                                 self._riskReversal25DeltaVols[i]*100.0)

            s += label_to_string("ATM Strike", self._K_ATM[i])
            s += label_to_string("ATM Delta", self._deltaATM[i])

            s += label_to_string("K_ATM", self._K_ATM[i])
            s += label_to_string("MS 25D Call Strike", self._K_25D_C_MS[i])
            s += label_to_string("MS 25D Put Strike", self._K_25D_P_MS[i])
            s += label_to_string("SKEW 25D CALL STRIKE", self._K_25D_C[i])
            s += label_to_string("SKEW 25D PUT STRIKE", self._K_25D_P[i])
            s += label_to_string("PARAMS", self._parameters[i])

        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
