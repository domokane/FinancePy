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
from ...utils.helpers import check_argument_types, label_to_string

from ...models.volatility_fns import VolFuncTypes
from ...models.volatility_fns import vol_function_clark
from ...models.volatility_fns import vol_function_bloomberg
from ...models.volatility_fns import vol_function_svi
from ...models.volatility_fns import vol_function_ssvi
from ...models.sabr import vol_function_sabr
from ...models.sabr import vol_function_sabr_beta_half
from ...models.sabr import vol_function_sabr_beta_one


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

# @njit(fastmath=True, cache=True)
# def _g(K, *args):
#     """ This is the objective function used in the determination of the FX
#     option implied strike which is computed in the class below. """

#     s = args[0]
#     t = args[1]
#     rd = args[2]
#     rf = args[3]
#     volatility = args[4]
#     delta_method_value = args[5]
#     option_type_value = args[6]
#     delta_target = args[7]

#     delta_out = fast_delta(s, t, K, rd, rf,
#                          volatility,
#                          delta_method_value,
#                          option_type_value)

#     obj_fn = delta_target - delta_out
#     return obj_fn

###############################################################################

# @njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
# def _interpolate_gap(k, strikes, gaps):

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
    """Return a value that is minimised when the ATM, MS and RR vols have
    been best fitted using the parametric volatility curve represented by
    params and specified by the vol_type_value at a single time slice only.
    """

    t = args[0]
    f = args[1]
    strikes_grid = args[2]
    index = args[3]
    vol_grid = args[4]
    vol_type_value = args[5]

    tot = 0.0

    num_strikes = len(vol_grid)

    for i in range(0, num_strikes):

        k = strikes_grid[i][index]
        fitted_vol = vol_function(vol_type_value, params, f, k, t)
        mkt_vol = vol_grid[i][index]
        diff = fitted_vol - mkt_vol
        tot += diff**2

    return tot


###############################################################################
# Do not cache this function as it leads to complaints
###############################################################################


def _solve_to_horizon(
    t,
    f,
    strikes_grid,
    time_index,
    vol_grid,
    vol_type_value,
    x_inits,
    fin_solver_type,
):

    ###########################################################################
    # Determine parameters of vol surface using minimisation
    ###########################################################################

    tol = 1e-6

    args = (t, f, strikes_grid, time_index, vol_grid, vol_type_value)

    # Nelder-Mead (both SciPy amd Numba) is quicker, but occasionally fails
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

    # print("t: %9.5f alpha:%9.5f beta: %9.5f rho: %9.5f nu: %9.5f" %
    #      (t, xopt[0], 0.5, xopt[1], xopt[2]))

    params = np.array(xopt)
    return params


###############################################################################


@njit(
    float64(int64, float64[:], float64, float64, float64),
    cache=True,
    fastmath=True,
)
def vol_function(vol_function_type_value, params, f, k, t):
    """Return the volatility for a strike using a given polynomial
    interpolation following Section 3.9 of Iain Clark book."""

    if vol_function_type_value == VolFuncTypes.CLARK.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.SABR_BETA_ONE.value:
        vol = vol_function_sabr_beta_one(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.SABR_BETA_HALF.value:
        vol = vol_function_sabr_beta_half(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.BBG.value:
        vol = vol_function_bloomberg(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.SABR.value:
        vol = vol_function_sabr(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.CLARK5.value:
        vol = vol_function_clark(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.SVI.value:
        vol = vol_function_svi(params, f, k, t)
        return vol
    elif vol_function_type_value == VolFuncTypes.SSVI.value:
        vol = vol_function_ssvi(params, f, k, t)
        return vol
    else:
        raise FinError("Unknown Model Type")


###############################################################################


# @njit(cache=True, fastmath=True)
# def _delta_fit(k, *args):
#     """ This is the objective function used in the determination of the FX
#     Option implied strike which is computed in the class below. I map it into
#     inverse normcdf space to avoid the flat slope of this function at low vol
#     and high K. It speeds up the code as it allows initial values close to
#     the solution to be used. """

#     vol_type_value = args[0]
#     s = args[1]
#     t = args[2]
#     rd = args[3]
#     rf = args[4]
#     option_type_value = args[5]
#     deltaTypeValue = args[6]
#     inverse_delta_target = args[7]
#     params = args[8]
#     strikes = args[9]
#     gaps = args[10]

#     f = s * np.exp((rd-rf)*t)
#     v = vol_function(vol_type_value, params, strikes, gaps, f, k, t)
#     delta_out = fast_delta(s, t, k, rd, rf, v, deltaTypeValue, option_type_value)
#     inverse_delta_out = norminvcdf(np.abs(delta_out))
#     inv_obj_fn = inverse_delta_target - inverse_delta_out

#     return inv_obj_fn

###############################################################################
# Unable to cache this function due to dynamic globals warning. Revisit.
###############################################################################


# @njit(float64(float64, float64, float64, float64, int64, int64, float64,
#              int64, float64, float64[:], float64[:], float64[:]),
#      fastmath=True)
# def _solver_for_smile_strike(s, t, rd, rf,
#                             option_type_value,
#                             volatilityTypeValue,
#                             delta_target,
#                             delta_method_value,
#                             initial_guess,
#                             parameters,
#                             strikes,
#                             gaps):
#     """ Solve for the strike that sets the delta of the option equal to the
#     target value of delta allowing the volatility to be a function of the
#     strike. """

#     inverse_delta_target = norminvcdf(np.abs(delta_target))

#     argtuple = (volatilityTypeValue, s, t, rd, rf,
#                 option_type_value, delta_method_value,
#                 inverse_delta_target,
#                 parameters, strikes, gaps)

#     K = newton_secant(_delta_fit, x0=initial_guess, args=argtuple,
#                       tol=1e-8, maxiter=50)

#     return K

###############################################################################
# Unable to cache function and if I remove njit it complains about pickle
###############################################################################


# @njit(float64(float64, float64, float64, float64, int64, float64,
#               int64, float64), fastmath=True)
# def solve_for_strike(spot_fx_rate,
#                    t_del, rd, rf,
#                    option_type_value,
#                    delta_target,
#                    delta_method_value,
#                    volatility):
#     """ This function determines the implied strike of an FX option
#     given a delta and the other option details. It uses a one-dimensional
#     Newton root search algorithm to determine the strike that matches an
#     input volatility. """

#     # ======================================================================
#     # IMPORTANT NOTE:
#     # ======================================================================
#     # For some delta quotation conventions I can solve for K explicitly.
#     # Note that as I am using the function norm_inv_delta to calculate the
#     # inverse value of delta, this may not, on a round trip using N(x), give
#     # back the value x as it is calculated to a different number of decimal
#     # places. It should however agree to 6-7 decimal places. Which is OK.
#     # ======================================================================

#     if delta_method_value == FinFXDeltaMethod.SPOT_DELTA.value:

#         dom_df = np.exp(-rd*t_del)
#         for_df = np.exp(-rf*t_del)

#         if option_type_value == OptionTypes.EUROPEAN_CALL.value:
#             phi = +1.0
#         else:
#             phi = -1.0

#         F0T = spot_fx_rate * for_df / dom_df
#         vsqrtt = volatility * np.sqrt(t_del)
#         arg = delta_target*phi/for_df  # CHECK THIS !!!
#         norm_inv_delta = norminvcdf(arg)
#         K = F0T * np.exp(-vsqrtt * (phi * norm_inv_delta - vsqrtt/2.0))
#         return K

#     elif delta_method_value == FinFXDeltaMethod.FORWARD_DELTA.value:

#         dom_df = np.exp(-rd*t_del)
#         for_df = np.exp(-rf*t_del)

#         if option_type_value == OptionTypes.EUROPEAN_CALL.value:
#             phi = +1.0
#         else:
#             phi = -1.0

#         F0T = spot_fx_rate * for_df / dom_df
#         vsqrtt = volatility * np.sqrt(t_del)
#         arg = delta_target*phi
#         norm_inv_delta = norminvcdf(arg)
#         K = F0T * np.exp(-vsqrtt * (phi * norm_inv_delta - vsqrtt/2.0))
#         return K

#     elif delta_method_value == FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ.value:

#         argtuple = (spot_fx_rate, t_del, rd, rf, volatility,
#                     delta_method_value, option_type_value, delta_target)

#         K = newton_secant(_g, x0=spot_fx_rate, args=argtuple,
#                           tol=1e-7, maxiter=50)

#         return K

#     elif delta_method_value == FinFXDeltaMethod.FORWARD_DELTA_PREM_ADJ.value:

#         argtuple = (spot_fx_rate, t_del, rd, rf, volatility,
#                     delta_method_value, option_type_value, delta_target)

#         K = newton_secant(_g, x0=spot_fx_rate, args=argtuple,
#                           tol=1e-7, maxiter=50)

#         return K

#     else:

#         raise FinError("Unknown FinFXDeltaMethod")

###############################################################################


class SwaptionVolSurface:
    """Class to perform a calibration of a chosen parametrised surface to the
    prices of swaptions at different expiry dates and swap tenors. There is a
    choice of volatility function from cubic in delta to full SABR and SSVI.
    Check out VolFuncTypes. Visualising the volatility curve is useful.
    Also, there is no guarantee that the implied pdf will be positive."""

    def __init__(
        self,
        value_dt: Date,
        expiry_dts: list,
        fwd_swap_rates: (list, np.ndarray),
        strike_grid: np.ndarray,
        vol_grid: np.ndarray,
        vol_func_type: VolFuncTypes = VolFuncTypes.SABR,
        fin_solver_type: FinSolverTypes = FinSolverTypes.NELDER_MEAD,
    ):
        """Create the FinSwaptionVolSurface object by passing in market vol
        data for a list of strikes and expiry dates."""

        check_argument_types(self.__init__, locals())

        self.value_dt = value_dt

        if len(strike_grid.shape) != 2:
            raise FinError("Strike grid must be a 2D grid of values")

        if len(vol_grid.shape) != 2:
            raise FinError("Volatility grid must be a 2D grid of values")

        if len(strike_grid) != len(vol_grid):
            raise FinError(
                "Strike grid and volatility grid must have same size"
            )

        if len(strike_grid[0]) != len(vol_grid[0]):
            raise FinError(
                "Strike grid and volatility grid must have same size"
            )

        if len(expiry_dts) != len(vol_grid[0]):
            raise FinError("Expiry dates not same size as volatility grid")

        self._num_expiry_dts = len(vol_grid[0])
        self._num_strikes = len(vol_grid)

        self._strike_grid = strike_grid
        self._vol_grid = vol_grid

        self._expiry_dts = expiry_dts
        self._vol_func_type = vol_func_type

        self._fwd_swap_rates = fwd_swap_rates

        self._build_vol_surface(fin_solver_type=fin_solver_type)

    #        self._F0T = []
    #        self._stock_price = None
    #        self._atm_method = None
    #        self._atm_vols = []
    #        self._delta_method = None
    #        self._strikes = []
    #        self._vol_grid = []

    ###########################################################################

    def vol_from_strike_dt(self, K, expiry_dt):
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

        vol_type_value = self._vol_func_type.value

        index0 = 0  # lower index in bracket
        index1 = 0  # upper index in bracket

        num_curves = self._num_expiry_dts

        if num_curves == 1:

            index0 = 0
            index1 = 0

        # If the time is below first time then assume a flat vol
        elif t_exp <= self._t_exp[0]:

            index0 = 0
            index1 = 0

        # If the time is beyond the last time then extrapolate with a flat vol
        elif t_exp >= self._t_exp[-1]:

            index0 = len(self._t_exp) - 1
            index1 = len(self._t_exp) - 1

        else:  # Otherwise we look for bracketing times and interpolate

            for i in range(1, num_curves):

                if t_exp <= self._t_exp[i] and t_exp > self._t_exp[i - 1]:
                    index0 = i - 1
                    index1 = i
                    break

        fwd0 = self._fwd_swap_rates[index0]
        fwd1 = self._fwd_swap_rates[index1]

        t0 = self._t_exp[index0]
        t1 = self._t_exp[index1]

        vol0 = vol_function(
            vol_type_value, self._parameters[index0], fwd0, K, t0
        )

        if index1 != index0:

            vol1 = vol_function(
                vol_type_value, self._parameters[index1], fwd1, K, t1
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

    # def delta_to_strike(self, call_delta, expiry_dt, delta_method):
    #     """ Interpolates the strike at a delta and expiry date. Linear
    #     interpolation is used in strike."""

    #     t_exp = (expiry_dt - self.value_dt) / g_days_in_year

    #     vol_type_value = self._vol_func_type.value

    #     s = self._spot_fx_rate

    #     if delta_method is None:
    #         delta_method_value = self._delta_method.value
    #     else:
    #         delta_method_value = delta_method.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     num_curves = self._num_vol_curves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if num_curves == 1:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is below first time then assume a flat vol
    #     elif t_exp <= self._t_exp[0]:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is beyond last time then extrapolate with a flat vol
    #     elif t_exp > self._t_exp[-1]:

    #         index0 = len(self._t_exp) - 1
    #         index1 = len(self._t_exp) - 1

    #     else: # Otherwise we look for bracketing times and interpolate

    #         for i in range(1, num_curves):

    #             if t_exp <= self._t_exp[i] and t_exp > self._t_exp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break

    #     #####################################################################

    #     t0 = self._t_exp[index0]
    #     t1 = self._t_exp[index1]

    #     initial_guess = self._k_atm[index0]

    #     K0 = _solver_for_smile_strike(s, t_exp, self._rd[index0],
    #                               self._rf[index0],
    #                               OptionTypes.EUROPEAN_CALL.value,
    #                               vol_type_value, call_delta,
    #                               delta_method_value,
    #                               initial_guess,
    #                               self._parameters[index0],
    #                               self._strikes[index0],
    #                               self._gaps[index0])

    #     if index1 != index0:

    #         K1 = _solver_for_smile_strike(s, t_exp,
    #                                   self._rd[index1],
    #                                   self._rf[index1],
    #                                   OptionTypes.EUROPEAN_CALL.value,
    #                                   vol_type_value, call_delta,
    #                                   delta_method_value,
    #                                   initial_guess,
    #                                   self._parameters[index1],
    #                                   self._strikes[index1],
    #                                   self._gaps[index1])
    #     else:

    #         K1 = K0

    #     # In the expiry time dimension, both volatilities are interpolated
    #     # at the same strikes but different deltas.

    #     if np.abs(t1-t0) > 1e-6:

    #         K = ((t_exp-t0) * K1 + (t1-t_exp) * K1) / (K1 - K0)

    #     else:

    #         K = K1

    #     return K

    ###########################################################################

    # def vol_from_delta_date(self, call_delta, expiry_dt,
    #                                delta_method = None):
    #     """ Interpolates the Black-Scholes volatility from the volatility
    #     surface given a call option delta and expiry date. Linear interpolation
    #     is done in variance space. The smile strikes at bracketed dates are
    #     determined by determining the strike that reproduces the provided delta
    #     value. This uses the calibration delta convention, but it can be
    #     overriden by a provided delta convention. The resulting volatilities
    #     are then determined for each bracketing expiry time and linear
    #     interpolation is done in variance space and then converted back to a
    #     lognormal volatility."""

    #     t_exp = (expiry_dt - self.value_dt) / g_days_in_year

    #     vol_type_value = self._vol_func_type.value

    #     s = self._spot_fx_rate

    #     if delta_method is None:
    #         delta_method_value = self._delta_method.value
    #     else:
    #         delta_method_value = delta_method.value

    #     index0 = 0 # lower index in bracket
    #     index1 = 0 # upper index in bracket

    #     num_curves = self._num_vol_curves

    #     # If there is only one time horizon then assume flat vol to this time
    #     if num_curves == 1:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is below first time then assume a flat vol
    #     elif t_exp <= self._t_exp[0]:

    #         index0 = 0
    #         index1 = 0

    #     # If the time is beyond the last time then extrapolate with a flat vol
    #     elif t_exp > self._t_exp[-1]:

    #         index0 = len(self._t_exp) - 1
    #         index1 = len(self._t_exp) - 1

    #     else: # Otherwise we look for bracketing times and interpolate

    #         for i in range(1, num_curves):

    #             if t_exp <= self._t_exp[i] and t_exp > self._t_exp[i-1]:
    #                 index0 = i-1
    #                 index1 = i
    #                 break

    #     fwd0 = self._F0T[index0]
    #     fwd1 = self._F0T[index1]

    #     t0 = self._t_exp[index0]
    #     t1 = self._t_exp[index1]

    #     initial_guess = self._k_atm[index0]

    #     K0 = _solver_for_smile_strike(s, t_exp, self._rd[index0], self._rf[index0],
    #                               OptionTypes.EUROPEAN_CALL.value,
    #                               vol_type_value, call_delta,
    #                               delta_method_value,
    #                               initial_guess,
    #                               self._parameters[index0],
    #                               self._strikes[index0],
    #                               self._gaps[index0])

    #     vol0 = vol_function(vol_type_value, self._parameters[index0],
    #                        self._strikes[index0], self._gaps[index0],
    #                        fwd0, K0, t0)

    #     if index1 != index0:

    #         K1 = _solver_for_smile_strike(s, t_exp,
    #                                   self._rd[index1],
    #                                   self._rf[index1],
    #                                   OptionTypes.EUROPEAN_CALL.value,
    #                                   vol_type_value, call_delta,
    #                                   delta_method_value,
    #                                   initial_guess,
    #                                   self._parameters[index1],
    #                                   self._strikes[index1],
    #                                   self._gaps[index1])

    #         vol1 = vol_function(vol_type_value, self._parameters[index1],
    #                            self._strikes[index1], self._gaps[index1],
    #                            fwd1, K1, t1)
    #     else:
    #         vol1 = vol0

    #     # In the expiry time dimension, both volatilities are interpolated
    #     # at the same strikes but different deltas.
    #     vart0 = vol0*vol0*t0
    #     vart1 = vol1*vol1*t1

    #     if np.abs(t1-t0) > 1e-6:

    #         vart = ((t_exp-t0) * vart1 + (t1-t_exp) * vart0) / (t1 - t0)
    #         kt = ((t_exp-t0) * K1 + (t1-t_exp) * K0) / (t1 - t0)

    #         if vart < 0.0:
    #             raise FinError("Failed interpolation due to negative variance.")

    #         volt = np.sqrt(vart/t_exp)

    #     else:

    #         volt = vol0
    #         kt = K0

    #     return volt, kt

    ###############################################################################

    def _build_vol_surface(self, fin_solver_type=FinSolverTypes.NELDER_MEAD):
        """Main function to construct the vol surface."""

        if self._vol_func_type == VolFuncTypes.CLARK:
            num_parameters = 3
        elif self._vol_func_type == VolFuncTypes.SABR_BETA_ONE:
            num_parameters = 3
        elif self._vol_func_type == VolFuncTypes.SABR_BETA_HALF:
            num_parameters = 3
        elif self._vol_func_type == VolFuncTypes.BBG:
            num_parameters = 3
        elif self._vol_func_type == VolFuncTypes.SABR:
            num_parameters = 4
        elif self._vol_func_type == VolFuncTypes.CLARK5:
            num_parameters = 5
        elif self._vol_func_type == VolFuncTypes.SVI:
            num_parameters = 5
        elif self._vol_func_type == VolFuncTypes.SSVI:
            num_parameters = 5
        else:
            print(self._vol_func_type)
            raise FinError("Unknown Model Type")

        num_expiry_dts = self._num_expiry_dts

        self._parameters = np.zeros([num_expiry_dts, num_parameters])
        self._t_exp = np.zeros(num_expiry_dts)

        #######################################################################
        # TODO: ADD SPOT DAYS
        #######################################################################

        for i in range(0, num_expiry_dts):

            expiry_dt = self._expiry_dts[i]
            t_exp = (expiry_dt - self.value_dt) / g_days_in_year
            self._t_exp[i] = t_exp

        #######################################################################
        # THE ACTUAL COMPUTATION LOOP STARTS HERE
        #######################################################################

        vol_type_value = self._vol_func_type.value

        x_inits = []
        x_init = np.zeros(num_parameters)
        x_inits.append(x_init)

        for i in range(0, num_expiry_dts):

            t = self._t_exp[i]
            f = self._fwd_swap_rates[i]

            res = _solve_to_horizon(
                t,
                f,
                self._strike_grid,
                i,
                self._vol_grid,
                vol_type_value,
                x_inits[i],
                fin_solver_type,
            )

            self._parameters[i, :] = res

            x_init = res
            x_inits.append(x_init)

    ###########################################################################

    def check_calibration(self, verbose: bool, tol: float = 1e-6):
        """Compare calibrated vol surface with market and output a report
        which sets out the quality of fit to the ATM and 10 and 25 delta market
        strangles and risk reversals."""

        if self._vol_grid == []:
            raise FinError("Error: Vol Grid is empty")

        if verbose:

            print("==========================================================")
            print("VALUE DATE:", self.value_dt)
            print("STOCK PRICE:", self._stock_price)
            print("==========================================================")

        for i in range(0, self._num_expiry_dts):

            expiry_dt = self._expiry_dts[i]
            print("==========================================================")

            for j in range(0, self._num_strikes):

                strike = self._strike_grid[j][i]
                fitted_vol = self.vol_from_strike_dt(strike, expiry_dt)
                mkt_vol = self._vol_grid[j][i]
                diff = fitted_vol - mkt_vol

                print(
                    "%s %12.3f %7.4f %7.4f %7.5f"
                    % (
                        expiry_dt,
                        strike,
                        fitted_vol * 100.0,
                        mkt_vol * 100,
                        diff * 100,
                    )
                )

        print("==========================================================")

    ###########################################################################

    # def implied_dbns(self, lowS, highS, num_intervals):
    #     """ Calculate the pdf for each tenor horizon. Returns a list of
    #     FinDistribution objects, one for each tenor horizon. """

    #     dbns = []

    #     for iTenor in range(0, self._num_expiry_dts):

    #         f = self._fwd_swap_rates[iTenor]
    #         t = self._t_exp[iTenor]

    #         dS = (highS - lowS)/ num_intervals

    #         dis_df = self._discount_curve.df(t)
    #         div_df = self._dividend_curve.df(t)

    #         r = -np.log(dis_df) / t
    #         q = -np.log(div_df) / t

    #         Ks = []
    #         vols = []

    #         for iK in range(0, num_intervals):

    #             k = lowS + iK*dS

    #             vol = vol_function(self._vol_func_type.value,
    #                               self._parameters[iTenor],
    #                               f, k, t)

    #             Ks.append(k)
    #             vols.append(vol)

    #         Ks = np.array(Ks)
    #         vols = np.array(vols)

    #         density = optionImpliedDbn(self._stock_price, t, r, q, Ks, vols)

    #         dbn = FinDistribution(Ks, density)
    #         dbns.append(dbn)

    #     return dbns

    ###########################################################################

    def plot_vol_curves(self):
        """Generates a plot of each of the vol discount implied by the market
        and fitted."""

        for tenor_index in range(0, self._num_expiry_dts):

            low_k = self._strike_grid[0][tenor_index] * 0.9
            high_k = self._strike_grid[-1][tenor_index] * 1.1

            expiry_dt = self._expiry_dts[tenor_index]
            plt.figure()

            ks = []

            num_intervals = 30
            K = low_k
            dK = (high_k - low_k) / num_intervals

            fitted_vols = []

            for i in range(0, num_intervals):

                ks.append(K)
                fitted_vol = self.vol_from_strike_dt(K, expiry_dt) * 100.0
                fitted_vols.append(fitted_vol)
                K = K + dK

            label_str = "FITTED AT " + str(self._expiry_dts[tenor_index])
            plt.plot(ks, fitted_vols, label=label_str)

            label_str = "MARKET AT " + str(self._expiry_dts[tenor_index])
            mkt_vols = self._vol_grid[:, tenor_index] * 100.0
            plt.plot(
                self._strike_grid[:, tenor_index],
                mkt_vols,
                "o",
                label=label_str,
            )

            plt.xlabel("Strike")
            plt.ylabel("Volatility")

            title = str(self._vol_func_type)
            plt.title(title)
            plt.legend()

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VALUE DATE", self.value_dt)
        s += label_to_string("STOCK PRICE", self._stock_price)
        s += label_to_string("ATM METHOD", self._atm_method)
        s += label_to_string("DELTA METHOD", self._delta_method)
        s += label_to_string("VOL FUNCTION", self._vol_func_type)

        for i in range(0, self._num_expiry_dts):

            s += "\n"

            s += label_to_string("EXPIRY DATE", self._expiry_dts[i])
            s += label_to_string("TIME (YRS)", self._t_exp[i])
            s += label_to_string("FWD FX", self._F0T[i])
            s += label_to_string("ATM VOLS", self._atm_vols[i] * 100.0)

            for j in range(0, self._num_strikes):

                expiry_dt = self._expiry_dts[i]
                k = self._strikes[j]
                vol = self._vol_grid[i][j]
                print(expiry_dt, k, vol)

        return s

    ###########################################################################

    def _print(self):
        """Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond."""
        print(self)


###############################################################################
