##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64, int64
from math import ceil

from ..utils.error import FinError
from ..utils.math import accrued_interpolator
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string
from ..utils.global_types import FinExerciseTypes
from ..utils.global_vars import gSmall

interp = InterpTypes.FLAT_FWD_RATES.value

###############################################################################
# TODO : Calculate accrued in bond option according to accrual convention
# TODO : Convergence is unstable - investigate how to improve it
# TODO : Write a fallback for gradient based alpha using bisection
# TODO : Fix treatment of accrued interest on the option expiry date.
###############################################################################

# (c) Dominic O'Kane - December-2019
# Fergal O'Kane - search root function - 16-12-2019

###############################################################################


def option_exercise_types_to_int(optionExerciseType):

    if optionExerciseType == FinExerciseTypes.EUROPEAN:
        return 1
    if optionExerciseType == FinExerciseTypes.BERMUDAN:
        return 2
    if optionExerciseType == FinExerciseTypes.AMERICAN:
        return 3
    else:
        raise FinError("Unknown option exercise type.")

###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def f(alpha, nm, Q, P, dX, dt, N):

    sumQZ = 0.0
    for j in range(-nm, nm+1):
        x = alpha + j*dX
        rdt = np.exp(x)*dt
        sumQZ += Q[j+N] * np.exp(-rdt)

    obj_fn = sumQZ - P
    return obj_fn

###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def fprime(alpha, nm, Q, P, dX, dt, N):

    sumQZdZ = 0.0
    for j in range(-nm, nm+1):
        x = alpha + j*dX
        rdt = np.exp(x)*dt
        sumQZdZ += Q[j+N] * np.exp(-rdt) * np.exp(x)

    deriv = -sumQZdZ*dt
    return deriv

###############################################################################
# This is the secant method which is not used as I computed the derivative of
# objective function with respect to the drift term
###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def search_root(x0, nm, Q, P, dX, dt, N):

    #    print("Searching for root", x0)
    max_iter = 50
    max_error = 1e-8

    x1 = x0 * 1.0001
    f0 = f(x0, nm, Q, P, dX, dt, N)
    f1 = f(x1, nm, Q, P, dX, dt, N)

    for _ in range(0, max_iter):

        df = f1 - f0

        if df == 0.0:
            raise FinError("Search for alpha fails due to zero derivative")

        x = x1 - f1 * (x1 - x0) / df
        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, nm, Q, P, dX, dt, N)

        if (abs(f1) <= max_error):
            return x1

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################
# This is Newton Raphson which is faster than the secant measure as it has the
# analytical derivative  that is easy to calculate.
###############################################################################


@njit(float64(float64, int64, float64[:], float64, float64, float64, int64),
      fastmath=True, cache=True)
def search_root_deriv(x0, nm, Q, P, dX, dt, N):

    max_iter = 50
    max_error = 1e-8

    for _ in range(0, max_iter):

        fval = f(x0, nm, Q, P, dX, dt, N)

        if abs(fval) <= max_error:
            return x0

        fderiv = fprime(x0, nm, Q, P, dX, dt, N)

        if abs(fderiv) == 0.0:
            print(x0, fval, fderiv)
            raise FinError("Function derivative is zero.")

        step = fval/fderiv
        x0 = x0 - step

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################


@njit(fastmath=True, cache=True)
def bermudan_swaption_tree_fast(texp, tmat,
                                strike_price, face_amount,
                                coupon_times, coupon_flows,
                                exercise_typeInt,
                                _df_times, _df_values,
                                _tree_times, _Q,
                                _pu, _pm, _pd,
                                _rt, _dt, _a):
    """ Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time. """

    num_time_steps, num_nodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################
    # I shove the floating rate value into the grid so it is handled in terms
    # of PV - it will not sit on a grid date so needs to be PV adjusted.
    # This is the value of the floating leg - this is a ONE CURVE approach.
    ###########################################################################

    fixed_legFlows = np.zeros(num_time_steps)
    # Initialise it with ones CHANGE
    float_leg_values = np.ones(num_time_steps)
    num_coupons = len(coupon_times)

    # swap fixed leg flows go all the way out to the swap maturity date
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        fixed_legFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree
        float_leg_values[n] = strike_price  # * df_flow / df_tree

    ############################## REMOVE START ###############################

    if 1 == 0:
        fixedpv = 0.0
        for n in range(0, num_time_steps-1):
            ttree = _tree_times[n]
            df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
            flow = fixed_legFlows[n]
            pvflow = flow * df_tree
            fixedpv += pvflow
            print(">>", n, ttree, df_tree, flow, fixedpv)
        fixedpv += df_tree
        df_tree = _uinterpolate(texp, _df_times, _df_values, interp)
        floatpv = df_tree
        swaptionpv = (fixedpv/df_tree - 1.0) * df_tree
        print("PV:", fixedpv, floatpv, swaptionpv)

        fixedpv = 0.0
        for n in range(0, num_coupons):
            tcpn = coupon_times[n]
            df = _uinterpolate(tcpn, _df_times, _df_values, interp)
            flow = coupon_flows[n]
            pvflow = flow * df
            fixedpv += pvflow
            print("++", n, tcpn, df, flow, fixedpv)
        fixedpv += df
        df_tree = _uinterpolate(texp, _df_times, _df_values, interp)
        floatpv = df_tree
        swaptionpv = (fixedpv/df_tree - 1.0) * df_tree
        print("PV:", fixedpv, floatpv, swaptionpv)

    ########################### REMOVE END ####################################

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mapped_times = np.array([0.0])
    mapped_amounts = np.array([0.0])
    for n in range(1, len(_tree_times)):

        accdAtExpiry = 0.0
        if _tree_times[n-1] < texp and _tree_times[n] >= texp:
            mapped_times = np.append(mapped_times, texp)
            mapped_amounts = np.append(mapped_amounts, accdAtExpiry)

        if fixed_legFlows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, fixed_legFlows[n])

    ###########################################################################

    accrued = np.zeros(num_time_steps)
    for m in range(0, maturityStep+1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, mapped_times, mapped_amounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if fixed_legFlows[m] > gSmall:
            accrued[m] = fixed_legFlows[m] * face_amount

    #######################################################################

    # The value of the swap at each time and node. Principal is exchanged.
    fixed_leg_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a payer swap
    pay_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a receiver swap
    rec_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Start with the value of the fixed leg at maturity
    for k in range(0, num_nodes):
        flow = 1.0 + fixed_legFlows[maturityStep]
        fixed_leg_values[maturityStep, k] = flow * face_amount

    N = jmax

    # Now step back to today considering early exercise on coupon dates
    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = fixed_legFlows[m] * face_amount

        for k in range(-nm, nm+1):
            kN = k + N
            rt = _rt[m, kN]
            df = np.exp(-rt * _dt)
            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = fixed_leg_values[m+1, kN]
                vm = fixed_leg_values[m+1, kN-1]
                vd = fixed_leg_values[m+1, kN-2]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixed_leg_values[m, kN] = v
            elif k == -jmax:
                vu = fixed_leg_values[m+1, kN+2]
                vm = fixed_leg_values[m+1, kN+1]
                vd = fixed_leg_values[m+1, kN]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixed_leg_values[m, kN] = v
            else:
                vu = fixed_leg_values[m+1, kN+1]
                vm = fixed_leg_values[m+1, kN]
                vd = fixed_leg_values[m+1, kN-1]
                v = (pu*vu + pm*vm + pd*vd) * df
                fixed_leg_values[m, kN] = v

            fixed_leg_values[m, kN] += flow
            vpay = 0.0
            vrec = 0.0

            if k == jmax:
                vu = pay_values[m+1, kN]
                vm = pay_values[m+1, kN-1]
                vd = pay_values[m+1, kN-2]
                vpay = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = pay_values[m+1, kN+2]
                vm = pay_values[m+1, kN+1]
                vd = pay_values[m+1, kN]
                vpay = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = pay_values[m+1, kN+1]
                vm = pay_values[m+1, kN]
                vd = pay_values[m+1, kN-1]
                vpay = (pu*vu + pm*vm + pd*vd) * df

            pay_values[m, kN] = vpay

            if k == jmax:
                vu = rec_values[m+1, kN]
                vm = rec_values[m+1, kN-1]
                vd = rec_values[m+1, kN-2]
                vrec = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = rec_values[m+1, kN+2]
                vm = rec_values[m+1, kN+1]
                vd = rec_values[m+1, kN]
                vrec = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = rec_values[m+1, kN+1]
                vm = rec_values[m+1, kN]
                vd = rec_values[m+1, kN-1]
                vrec = (pu*vu + pm*vm + pd*vd) * df

            rec_values[m, kN] = vrec

            holdPay = pay_values[m, kN]
            holdRec = rec_values[m, kN]

            # The floating value is clean and so must be the fixed value
            fixed_leg_value = fixed_leg_values[m, kN] - accrued[m]
            float_leg_value = float_leg_values[m]

            payExercise = max(float_leg_value - fixed_leg_value, 0.0)
            recExercise = max(fixed_leg_value - float_leg_value, 0.0)

            if m == expiryStep:

                pay_values[m, kN] = max(payExercise, holdPay)
                rec_values[m, kN] = max(recExercise, holdRec)

            elif exercise_typeInt == 2 and flow > gSmall and m > expiryStep:

                pay_values[m, kN] = max(payExercise, holdPay)
                rec_values[m, kN] = max(recExercise, holdRec)

            elif exercise_typeInt == 3 and m > expiryStep:

                raise FinError("American optionality not allowed.")

                # Need to define floating value on all grid dates

                pay_values[m, kN] = max(payExercise, holdPay)
                rec_values[m, kN] = max(recExercise, holdRec)

    return pay_values[0, jmax], rec_values[0, jmax]

###############################################################################


@njit(fastmath=True, cache=True)
def american_bond_option_tree_fast(texp, tmat,
                                   strike_price, face_amount,
                                   coupon_times, coupon_flows,
                                   exercise_typeInt,
                                   _df_times, _df_values,
                                   _tree_times, _Q,
                                   _pu, _pm, _pd,
                                   _rt,
                                   _dt, _a):
    """ Option that can be exercised at any time over the exercise period.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. """

    DEBUG = False

    ###########################################################################

    num_time_steps, num_nodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    treeFlows = np.zeros(num_time_steps)
    num_coupons = len(coupon_times)

    # Tree flows go all the way out to the bond maturity date
    # Do not include first coupon as it is the previous coupon and is negative
    for i in range(1, num_coupons):
        tcpn = coupon_times[i]

        if tcpn < 0.0:
            raise FinError("Coupon times must be positive.")

        n = int(tcpn/_dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        treeFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    #
    # mapped_times = np.zeros(0)
    # mapped_amounts = np.zeros(0)
    # for n in range(0, len(_tree_times)):
    #     if treeFlows[n] > 0.0:
    #         mapped_times = np.append(mapped_times, _tree_times[n])
    #         mapped_amounts = np.append(mapped_amounts, treeFlows[n])
    ###########################################################################

    accrued = np.zeros(num_time_steps)
    for m in range(0, maturityStep+1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, coupon_times, coupon_flows)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > gSmall:
            accrued[m] = treeFlows[m] * face_amount

    if DEBUG:
        for i in range(0, expiryStep+1):
            print(i, treeFlows[i], accrued[i])

    #######################################################################

    call_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    put_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Start with the value of the bond at maturity
    for k in range(0, num_nodes):
        bond_values[maturityStep, k] = (1.0 + treeFlows[maturityStep]) \
            * face_amount

    if DEBUG:
        full_price = bond_values[maturityStep, 0]
        clean_price = full_price - accrued[maturityStep]
        print(m, _tree_times[m], accrued[m], full_price, clean_price, 0, 0)

    # Step back from maturity to expiry date but with no exercise allowed.
    for m in range(maturityStep-1, expiryStep, -1):

        nm = min(m, jmax)
        flow = treeFlows[m] * face_amount

        for k in range(-nm, nm+1):
            kN = k + jmax
            r = _rt[m, kN]
            df = np.exp(-r * _dt)

            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = bond_values[m+1, kN]
                vm = bond_values[m+1, kN-1]
                vd = bond_values[m+1, kN-2]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v
            elif k == -jmax:
                vu = bond_values[m+1, kN+2]
                vm = bond_values[m+1, kN+1]
                vd = bond_values[m+1, kN]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v
            else:
                vu = bond_values[m+1, kN+1]
                vm = bond_values[m+1, kN]
                vd = bond_values[m+1, kN-1]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v

            bond_values[m, kN] += flow

        if DEBUG:
            print(m, _tree_times[m], accrued[m], full_price, clean_price, 0, 0)

    # Now consider exercise of the option on and before the expiry date
    for m in range(expiryStep, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face_amount

        for k in range(-nm, nm+1):

            kN = k + jmax
            r = _rt[m, kN]
            df = np.exp(-r * _dt)

            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = bond_values[m+1, kN]
                vm = bond_values[m+1, kN-1]
                vd = bond_values[m+1, kN-2]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v
            elif k == -jmax:
                vu = bond_values[m+1, kN+2]
                vm = bond_values[m+1, kN+1]
                vd = bond_values[m+1, kN]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v
            else:
                vu = bond_values[m+1, kN+1]
                vm = bond_values[m+1, kN]
                vd = bond_values[m+1, kN-1]
                v = (pu*vu + pm*vm + pd*vd) * df
                bond_values[m, kN] = v

            bond_values[m, kN] += flow

            vcall = 0.0
            vput = 0.0

            if k == jmax:
                vu = call_option_values[m+1, kN]
                vm = call_option_values[m+1, kN-1]
                vd = call_option_values[m+1, kN-2]
                vcall = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = call_option_values[m+1, kN+2]
                vm = call_option_values[m+1, kN+1]
                vd = call_option_values[m+1, kN]
                vcall = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = call_option_values[m+1, kN+1]
                vm = call_option_values[m+1, kN]
                vd = call_option_values[m+1, kN-1]
                vcall = (pu*vu + pm*vm + pd*vd) * df

            call_option_values[m, kN] = vcall

            if k == jmax:
                vu = put_option_values[m+1, kN]
                vm = put_option_values[m+1, kN-1]
                vd = put_option_values[m+1, kN-2]
                vput = (pu*vu + pm*vm + pd*vd) * df
            elif k == -jmax:
                vu = put_option_values[m+1, kN+2]
                vm = put_option_values[m+1, kN+1]
                vd = put_option_values[m+1, kN]
                vput = (pu*vu + pm*vm + pd*vd) * df
            else:
                vu = put_option_values[m+1, kN+1]
                vm = put_option_values[m+1, kN]
                vd = put_option_values[m+1, kN-1]
                vput = (pu*vu + pm*vm + pd*vd) * df

            put_option_values[m, kN] = vput

            full_price = bond_values[m, kN]
            clean_price = full_price - accrued[m]
            callExercise = max(clean_price - strike_price, 0.0)
            putExercise = max(strike_price - clean_price, 0.0)

            holdCall = call_option_values[m, kN]
            holdPut = put_option_values[m, kN]

            if m == expiryStep:

                call_option_values[m, kN] = max(callExercise, holdCall)
                put_option_values[m, kN] = max(putExercise, holdPut)

            elif exercise_typeInt == 3 and m < expiryStep:  # AMERICAN

                call_option_values[m, kN] = max(callExercise, holdCall)
                put_option_values[m, kN] = max(putExercise, holdPut)

        if DEBUG:
            print(m, _tree_times[m], accrued[m], full_price, clean_price,
                  callExercise, putExercise)

    return call_option_values[0, jmax], put_option_values[0, jmax]

###############################################################################


@njit(fastmath=True, cache=True)
def callable_puttable_bond_tree_fast(coupon_times, coupon_flows,
                                     call_times, call_prices,
                                     put_times, put_prices, face_amount,
                                     _sigma, _a, _Q,  # IS SIGMA USED ?
                                     _pu, _pm, _pd, _rt, _dt, _tree_times,
                                     _df_times, _df_values):
    """ Value a bond with embedded put and call options that can be exercised
    at any time over the specified list of put and call dates.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. """

    #######################################################################
    num_time_steps, num_nodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * _dt))
    tmat = coupon_times[-1]
    maturityStep = int(tmat/dt + 0.50)

    ###########################################################################
    # Map coupons onto tree while preserving their present value
    ###########################################################################

    treeFlows = np.zeros(num_time_steps)

    num_coupons = len(coupon_times)
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        treeFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree

    #######################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    #######################################################################

    mapped_times = np.array([0.0])
    mapped_amounts = np.array([0.0])
    for n in range(1, len(_tree_times)):
        if treeFlows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, treeFlows[n])

    #######################################################################

    accrued = np.zeros(num_time_steps)
    for m in range(0, num_time_steps):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, mapped_times, mapped_amounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face_amount

    ###########################################################################
    # map call onto tree - must have no calls at high value
    ###########################################################################

    tree_call_value = np.ones(num_time_steps) * face_amount * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(call_time/dt + 0.50)
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    treePutValue = np.zeros(num_time_steps)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(put_time/dt + 0.50)
        treePutValue[n] = put_prices[i]

    ###########################################################################
    # Value the bond by backward induction starting at bond maturity
    ###########################################################################

    callPutBondValues = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    DEBUG = False
    if DEBUG:
        df = 1.0
        px = 0.0
        for i in range(0, maturityStep+1):
            flow = treeFlows[i]
            t = _tree_times[i]
            df = _uinterpolate(t, _df_times, _df_values, interp)
            px += flow*df
        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = min(maturityStep, jmax)
    vcall = tree_call_value[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * face_amount
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm+1):
        kN = k + jmax
        bond_values[m, kN] = (1.0 + treeFlows[m]) * face_amount
        callPutBondValues[m, kN] = value

    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face_amount
        vcall = tree_call_value[m]
        vput = treePutValue[m]

        for k in range(-nm, nm+1):
            kN = k + jmax
            rt = _rt[m, kN]
            df = np.exp(-rt*dt)
            pu = _pu[kN]
            pm = _pm[kN]
            pd = _pd[kN]

            if k == jmax:
                vu = bond_values[m+1, kN]
                vm = bond_values[m+1, kN-1]
                vd = bond_values[m+1, kN-2]
            elif k == -jmax:
                vu = bond_values[m+1, kN+2]
                vm = bond_values[m+1, kN+1]
                vd = bond_values[m+1, kN]
            else:
                vu = bond_values[m+1, kN+1]
                vm = bond_values[m+1, kN]
                vd = bond_values[m+1, kN-1]

            v = (pu*vu + pm*vm + pd*vd) * df
            bond_values[m, kN] = v
            bond_values[m, kN] += flow

            if k == jmax:
                vu = callPutBondValues[m+1, kN]
                vm = callPutBondValues[m+1, kN-1]
                vd = callPutBondValues[m+1, kN-2]
            elif k == -jmax:
                vu = callPutBondValues[m+1, kN+2]
                vm = callPutBondValues[m+1, kN+1]
                vd = callPutBondValues[m+1, kN]
            else:
                vu = callPutBondValues[m+1, kN+1]
                vm = callPutBondValues[m+1, kN]
                vd = callPutBondValues[m+1, kN-1]

            vhold = (pu*vu + pm*vm + pd*vd) * df
            # Need to make add on coupons paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            callPutBondValues[m, kN] = value

    return {'bondwithoption': callPutBondValues[0, jmax],
            'bondpure': bond_values[0, jmax]}

###############################################################################


@njit(fastmath=True, cache=True)
def build_tree_fast(a, sigma, tree_times, num_time_steps, discount_factors):
    """ Calibrate the tree to a term structure of interest rates. """

    treeMaturity = tree_times[-1]
    dt = treeMaturity / (num_time_steps+1)
    dX = sigma * np.sqrt(3.0 * dt)
    jmax = ceil(0.1835/(a * dt))

    if jmax > 1000:
        raise FinError("Jmax > 1000. Increase a or dt.")

    pu = np.zeros(shape=(2*jmax+1))
    pm = np.zeros(shape=(2*jmax+1))
    pd = np.zeros(shape=(2*jmax+1))

    # The short rate goes out one step extra to have the final short rate
    # This is the BK model so x = log(r)
    X = np.zeros(shape=(num_time_steps+2, 2*jmax+1))
    rt = np.zeros(shape=(num_time_steps+2, 2*jmax+1))

    # probabilities start at time 0 and go out to one step before T
    # Branching is simple trinomial out to time step m=1 after which
    # the top node and bottom node connect internally to two lower nodes
    # and two upper nodes respectively. The probabilities only depend on j

    for j in range(-jmax, jmax+1):

        ajdt = a*j*dt
        jN = j + jmax

        if j == jmax:
            pu[jN] = 7.0/6.0 + 0.50 * (ajdt*ajdt - 3.0*ajdt)
            pm[jN] = -1.0/3.0 - ajdt * ajdt + 2.0 * ajdt
            pd[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt - ajdt)
        elif j == -jmax:
            pu[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt + ajdt)
            pm[jN] = -1.0/3.0 - ajdt * ajdt - 2.0 * ajdt
            pd[jN] = 7.0/6.0 + 0.50 * (ajdt * ajdt + 3.0 * ajdt)
        else:
            pu[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt - ajdt)
            pm[jN] = 2.0/3.0 - ajdt * ajdt
            pd[jN] = 1.0/6.0 + 0.50 * (ajdt * ajdt + ajdt)

    # Arrow-Debreu array
    Q = np.zeros(shape=(num_time_steps+2, 2*jmax+1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(num_time_steps+1)

    # Time zero is trivial for the Arrow-Debreu price
    Q[0, jmax] = 1.0

    # Estimate short rate over first year
    r0 = -np.log(discount_factors[1])/tree_times[1]

    # We initialise x0 with value of log of r0
    x0 = np.log(r0)

    # Big loop over time steps
    for m in range(0, num_time_steps + 1):

        nm = min(m, jmax)

        # Need to do drift adjustment which is non-linear and so requires
        # a root search algorithm to find value of x0.

        alpha[m] = search_root_deriv(x0, nm, Q[m], discount_factors[m + 1],
                                     dX, dt, jmax)

        x0 = alpha[m]

        for j in range(-nm, nm+1):
            jN = j + jmax
            X[m, jN] = alpha[m] + j*dX
            rt[m, jN] = np.exp(X[m, jN])

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):

            jN = j + jmax
            rdt = np.exp(X[m, jN]) * dt
            z = np.exp(-rdt)

            if j == jmax:
                Q[m+1, jN] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-2] += Q[m, jN] * pd[jN] * z
            elif j == -jmax:
                Q[m+1, jN] += Q[m, jN] * pd[jN] * z
                Q[m+1, jN+1] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN+2] += Q[m, jN] * pu[jN] * z
            else:
                Q[m+1, jN+1] += Q[m, jN] * pu[jN] * z
                Q[m+1, jN] += Q[m, jN] * pm[jN] * z
                Q[m+1, jN-1] += Q[m, jN] * pd[jN] * z

    return (Q, pu, pm, pd, rt, dt)

##########################################################################


class BKTree():

    def __init__(self,
                 sigma: float,
                 a: float,
                 num_time_steps: int = 100):
        """ Constructs the Black Karasinski rate model. The speed of mean
        reversion a and volatility are passed in. The short rate process
        is given by d(log(r)) = (theta(t) - a*log(r)) * dt  + sigma * dW """

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        if a < 1e-10:
            a = 1e-10

        self._a = a
        self._sigma = sigma

        if num_time_steps < 3:
            raise FinError("Drift fitting requires at least 3 time steps")

        self._num_time_steps = num_time_steps

        self._Q = None
        self._rt = None
        self._tree_times = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discount_curve = None

###############################################################################

    def build_tree(self, tmat, df_times, df_values):

        if isinstance(df_times, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(df_values, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = InterpTypes.FLAT_FWD_RATES.value

        treeMaturity = tmat * (self._num_time_steps+1)/self._num_time_steps
        tree_times = np.linspace(0.0, treeMaturity, self._num_time_steps + 2)
        self._tree_times = tree_times

        dfTree = np.zeros(shape=(self._num_time_steps+2))
        dfTree[0] = 1.0

        for i in range(1, self._num_time_steps+2):
            t = tree_times[i]
            dfTree[i] = _uinterpolate(t, df_times, df_values, interp)

        self._df_times = df_times
        self._dfs = df_values

        self._Q, self._pu, self._pm, self._pd, self._rt, self._dt \
            = build_tree_fast(self._a, self._sigma,
                              tree_times, self._num_time_steps, dfTree)

        return

###############################################################################

    def bond_option(self, texp, strike_price, face_amount,
                    coupon_times, coupon_flows, exercise_type):
        """ Value a bond option that has European or American exercise using
        the Black-Karasinski model. The model uses a trinomial tree. """

        exercise_typeInt = option_exercise_types_to_int(exercise_type)

        tmat = coupon_times[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        callValue, putValue \
            = american_bond_option_tree_fast(texp, tmat,
                                             strike_price, face_amount,
                                             coupon_times, coupon_flows,
                                             exercise_typeInt,
                                             self._df_times, self._dfs,
                                             self._tree_times, self._Q,
                                             self._pu, self._pm, self._pd,
                                             self._rt,
                                             self._dt, self._a)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudan_swaption(self, texp, tmat, strike_price, face_amount,
                          coupon_times, coupon_flows, exercise_type):
        """ Swaption that can be exercised on specific dates over the exercise
        period. Due to non-analytical bond price we need to extend tree out to
        bond maturity and take into account cash flows through time. """

        exercise_typeInt = option_exercise_types_to_int(exercise_type)

        tmat = coupon_times[-1]

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        #######################################################################

        payValue, recValue \
            = bermudan_swaption_tree_fast(texp, tmat,
                                          strike_price, face_amount,
                                          coupon_times, coupon_flows,
                                          exercise_typeInt,
                                          self._df_times, self._dfs,
                                          self._tree_times, self._Q,
                                          self._pu, self._pm, self._pd,
                                          self._rt,
                                          self._dt, self._a)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def callable_puttable_bond_tree(self,
                                    coupon_times, coupon_flows,
                                    call_times, call_prices,
                                    put_times, put_prices,
                                    face):
        """ Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. """

        call_times = np.array(call_times)
        put_times = np.array(put_times)

        call_prices = np.array(call_prices)
        put_prices = np.array(put_prices)

        v = callable_puttable_bond_tree_fast(coupon_times, coupon_flows,
                                             call_times, call_prices,
                                             put_times, put_prices, face,
                                             self._sigma, self._a,
                                             self._Q,
                                             self._pu, self._pm, self._pd,
                                             self._rt, self._dt,
                                             self._tree_times,
                                             self._df_times, self._dfs)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def __repr__(self):
        """ Return string with class details. """

        s = "Black-Karasinski Model\n"
        s += label_to_string("Sigma", self._sigma)
        s += label_to_string("a", self._a)
        s += label_to_string("num_time_steps", self._num_time_steps)
        return s

###############################################################################
