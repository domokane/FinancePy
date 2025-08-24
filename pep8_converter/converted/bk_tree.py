# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from math import ceil

import numpy as np
from numba import njit, float64, int64

from ..utils.error import FinError
from ..utils.math import accrued_interpolator
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string
from ..utils.global_types import FinExerciseTypes
from ..utils.global_vars import G_SMALL

interp = InterpTypes.FLAT_FWD_RATES.value

# TODO : Calculate accrued in bond option according to accrual convention
# TODO : Convergence is unstable - investigate how to improve it
# TODO : Write a fallback for gradient based alpha using bisection
# TODO : Fix treatment of accrued interest on the option expiry date.

# (c) Dominic O'Kane - December-2019
# Fergal O'Kane - search root function - 16-12-2019

########################################################################################

def option_exercise_types_to_int(option_exercise_type):

    if option_exercise_type == FinExerciseTypes.EUROPEAN:
        return 1
    if option_exercise_type == FinExerciseTypes.BERMUDAN:
        return 2
    if option_exercise_type == FinExerciseTypes.AMERICAN:
        return 3
    else:
        raise FinError("Unknown option exercise type.")

########################################################################################

@njit(
    float64(float64, int64, float64[:], float64, float64, float64, int64),
    fastmath=True,
    cache=True,
)
def f(alpha, nm, qq, pp, dx, dt, n):

    sum_qz = 0.0
    for j in range(-nm, nm + 1):
        x = alpha + j * dx
        rdt = np.exp(x) * dt
        sum_qz += qq[j + n] * np.exp(-rdt)

    obj_fn = sum_qz - pp
    return obj_fn

########################################################################################

@njit(
    float64(float64, int64, float64[:], float64, float64, float64, int64),
    fastmath=True,
    cache=True,
)
def fprime(alpha, nm, qq, pp, dx, dt, n):

    """Compute the derivative of the objective function."""
    sum_q_zd_z = 0.0
    for j in range(-nm, nm + 1):
        x = alpha + j * dx
        rdt = np.exp(x) * dt
        sum_q_zd_z += qq[j + n] * np.exp(-rdt) * np.exp(x)

    deriv = -sum_q_zd_z * dt
    return deriv


# This is the secant method which is not used as I computed the derivative of
# objective function with respect to the drift term

########################################################################################

@njit(
    float64(float64, int64, float64[:], float64, float64, float64, int64),
    fastmath=True,
    cache=True,
)
def search_root(x0, nm, qq, pp, dx, dt, n):

    #    print("Searching for root", x0)
    max_iter = 50
    max_error = 1e-8

    x1 = x0 * 1.0001
    f0 = f(x0, nm, qq, pp, dx, dt, n)
    f1 = f(x1, nm, qq, pp, dx, dt, n)

    for _ in range(0, max_iter):

        df = f1 - f0

        if df == 0.0:
            raise FinError("Search for alpha fails due to zero derivative")

        x = x1 - f1 * (x1 - x0) / df
        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, nm, qq, pp, dx, dt, n)

        if abs(f1) <= max_error:
            return x1

    raise FinError("Search root deriv FAILED to find alpha.")


# This is Newton Raphson which is faster than the secant measure as it has the
# analytical derivative  that is easy to calculate.

########################################################################################

@njit(
    float64(float64, int64, float64[:], float64, float64, float64, int64),
    fastmath=True,
    cache=True,
)
def search_root_deriv(x0, nm, qq, pp, dx, dt, n):

    max_iter = 50
    max_error = 1e-8

    for _ in range(0, max_iter):

        fval = f(x0, nm, qq, pp, dx, dt, n)

        if abs(fval) <= max_error:
            return x0

        fderiv = fprime(x0, nm, qq, pp, dx, dt, n)

        if abs(fderiv) == 0.0:
            print(x0, fval, fderiv)
            raise FinError("Function derivative is zero.")

        step = fval / fderiv
        x0 = x0 - step

    raise FinError("Search root deriv FAILED to find alpha.")

########################################################################################

@njit(fastmath=True, cache=True)
def bermudan_swaption_tree_fast(

    t_exp,
    t_mat,
    strike_price,
    face_amount,
    cpn_times,
    cpn_flows,
    exercise_type_int,
    _df_times,
    _df_values,
    _tree_times,
    _qq,
    _pu,
    _pm,
    _pd,
    _rt,
    _dt,
    _a,
):
    """Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time."""

    num_time_steps, num_nodes = _qq.shape
    j_max = ceil(0.1835 / (_a * _dt))
    expiry_step = int(t_exp / _dt + 0.50)
    maturity_step = int(t_mat / _dt + 0.50)

    # I shove the floating rate value into the grid, so it is handled in terms
    # of ppV - it will not sit on a grid date so needs to be ppV adjusted.
    # This is the value of the floating leg - this is a ONE CURVE approach.

    fixed_leg_flows = np.zeros(num_time_steps)
    # Initialise it with ones CHANGE
    float_leg_values = np.ones(num_time_steps)
    num_cpns = len(cpn_times)

    # swap fixed leg flows go all the way out to the swap maturity date
    for i in range(0, num_cpns):
        t_cpn = cpn_times[i]
        n = int(t_cpn / _dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(t_cpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        fixed_leg_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree
        float_leg_values[n] = strike_price  # * df_flow / df_tree

    ############################## REMOVE START ###############################

    if 1 == 0:
        fixed_pv = 0.0
        for n in range(0, num_time_steps - 1):
            ttree = _tree_times[n]
            df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
            flow = fixed_leg_flows[n]
            pv_flow = flow * df_tree
            fixed_pv += pv_flow
            print(">>", n, ttree, df_tree, flow, fixed_pv)
        fixed_pv += df_tree
        df_tree = _uinterpolate(t_exp, _df_times, _df_values, interp)
        floatpv = df_tree
        swaptionpv = (fixed_pv / df_tree - 1.0) * df_tree
        print("ppV:", fixed_pv, floatpv, swaptionpv)

        fixed_pv = 0.0
        for n in range(0, num_cpns):
            t_cpn = cpn_times[n]
            df = _uinterpolate(t_cpn, _df_times, _df_values, interp)
            flow = cpn_flows[n]
            pv_flow = flow * df
            fixed_pv += pv_flow
            print("++", n, t_cpn, df, flow, fixed_pv)
        fixed_pv += df
        df_tree = _uinterpolate(t_exp, _df_times, _df_values, interp)
        floatpv = df_tree
        swaptionpv = (fixed_pv / df_tree - 1.0) * df_tree
        print("ppV:", fixed_pv, floatpv, swaptionpv)

    ########################### REMOVE END ####################################

    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent

    mapped_times = np.array([0.0])
    mapped_amounts = np.array([0.0])
    for n in range(1, len(_tree_times)):

        accd_at_expiry = 0.0
        if _tree_times[n - 1] < t_exp and _tree_times[n] >= t_exp:
            mapped_times = np.append(mapped_times, t_exp)
            mapped_amounts = np.append(mapped_amounts, accd_at_expiry)

        if fixed_leg_flows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, fixed_leg_flows[n])


    accrued = np.zeros(num_time_steps)
    for m in range(0, maturity_step + 1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, mapped_times, mapped_amounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if fixed_leg_flows[m] > G_SMALL:
            accrued[m] = fixed_leg_flows[m] * face_amount


    # The value of the swap at each time and node. pprincipal is exchanged.
    fixed_leg_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a payer swap
    pay_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a receiver swap
    rec_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Start with the value of the fixed leg at maturity
    for k in range(0, num_nodes):
        flow = 1.0 + fixed_leg_flows[maturity_step]
        fixed_leg_values[maturity_step, k] = flow * face_amount

    n = j_max

    # Now step back to today considering early exercise on coupon dates
    for m in range(maturity_step - 1, -1, -1):
        nm = min(m, j_max)
        flow = fixed_leg_flows[m] * face_amount

        for k in range(-nm, nm + 1):
            kn = k + n
            rt = _rt[m, kn]
            df = np.exp(-rt * _dt)
            pu = _pu[kn]
            pm = _pm[kn]
            pd = _pd[kn]

            if k == j_max:
                vu = fixed_leg_values[m + 1, kn]
                vm = fixed_leg_values[m + 1, kn - 1]
                vd = fixed_leg_values[m + 1, kn - 2]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, kn] = v
            elif k == -j_max:
                vu = fixed_leg_values[m + 1, kn + 2]
                vm = fixed_leg_values[m + 1, kn + 1]
                vd = fixed_leg_values[m + 1, kn]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, kn] = v
            else:
                vu = fixed_leg_values[m + 1, kn + 1]
                vm = fixed_leg_values[m + 1, kn]
                vd = fixed_leg_values[m + 1, kn - 1]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, kn] = v

            fixed_leg_values[m, kn] += flow
            vpay = 0.0
            vrec = 0.0

            if k == j_max:
                vu = pay_values[m + 1, kn]
                vm = pay_values[m + 1, kn - 1]
                vd = pay_values[m + 1, kn - 2]
                vpay = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = pay_values[m + 1, kn + 2]
                vm = pay_values[m + 1, kn + 1]
                vd = pay_values[m + 1, kn]
                vpay = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = pay_values[m + 1, kn + 1]
                vm = pay_values[m + 1, kn]
                vd = pay_values[m + 1, kn - 1]
                vpay = (pu * vu + pm * vm + pd * vd) * df

            pay_values[m, kn] = vpay

            if k == j_max:
                vu = rec_values[m + 1, kn]
                vm = rec_values[m + 1, kn - 1]
                vd = rec_values[m + 1, kn - 2]
                vrec = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = rec_values[m + 1, kn + 2]
                vm = rec_values[m + 1, kn + 1]
                vd = rec_values[m + 1, kn]
                vrec = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = rec_values[m + 1, kn + 1]
                vm = rec_values[m + 1, kn]
                vd = rec_values[m + 1, kn - 1]
                vrec = (pu * vu + pm * vm + pd * vd) * df

            rec_values[m, kn] = vrec

            hold_pay = pay_values[m, kn]
            hold_rec = rec_values[m, kn]

            # The floating value is clean and so must be the fixed value
            fixed_leg_value = fixed_leg_values[m, kn] - accrued[m]
            float_leg_value = float_leg_values[m]

            pay_exercise = max(float_leg_value - fixed_leg_value, 0.0)
            rec_exercise = max(fixed_leg_value - float_leg_value, 0.0)

            if m == expiry_step:

                pay_values[m, kn] = max(pay_exercise, hold_pay)
                rec_values[m, kn] = max(rec_exercise, hold_rec)

            elif exercise_type_int == 2 and flow > G_SMALL and m > expiry_step:

                pay_values[m, kn] = max(pay_exercise, hold_pay)
                rec_values[m, kn] = max(rec_exercise, hold_rec)

            elif exercise_type_int == 3 and m > expiry_step:

                raise FinError("American optionality not allowed.")

                # Need to define floating value on all grid dates

                # pay_values[m, kn] = max(pay_exercise, hold_pay)
                # rec_values[m, kn] = max(rec_exercise, hold_rec)

    return pay_values[0, j_max], rec_values[0, j_max]

########################################################################################

@njit(fastmath=True, cache=True)
def american_bond_option_tree_fast(

    t_exp,
    t_mat,
    strike_price,
    face_amount,
    cpn_times,
    cpn_flows,
    exercise_type_int,
    _df_times,
    _df_values,
    _tree_times,
    _qq,
    _pu,
    _pm,
    _pd,
    _rt,
    _dt,
    _a,
):
    """Option that can be exercised at any time over the exercise period.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time."""

    debug = False


    num_time_steps, num_nodes = _qq.shape
    j_max = ceil(0.1835 / (_a * _dt))
    expiry_step = int(t_exp / _dt + 0.50)
    maturity_step = int(t_mat / _dt + 0.50)


    tree_flows = np.zeros(num_time_steps)
    num_cpns = len(cpn_times)

    # Tree flows go all the way out to the bond maturity date
    # Do not include first coupon as it is the previous coupon and is negative
    for i in range(1, num_cpns):
        t_cpn = cpn_times[i]

        if t_cpn < 0.0:
            raise FinError("Coupon times must be positive.")

        n = int(t_cpn / _dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(t_cpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        tree_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree

    # mapped_times = np.zeros(0)
    # mapped_amounts = np.zeros(0)
    # for n in range(0, len(_tree_times)):
    #     if tree_flows[n] > 0.0:
    #         mapped_times = np.append(mapped_times, _tree_times[n])
    #         mapped_amounts = np.append(mapped_amounts, tree_flows[n])

    accrued = np.zeros(num_time_steps)
    for m in range(0, maturity_step + 1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, cpn_times, cpn_flows)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if tree_flows[m] > G_SMALL:
            accrued[m] = tree_flows[m] * face_amount

    if debug:
        for i in range(0, expiry_step + 1):
            print(i, tree_flows[i], accrued[i])


    call_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    put_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Start with the value of the bond at maturity
    for k in range(0, num_nodes):
        bond_values[maturity_step, k] = (1.0 + tree_flows[maturity_step]) * face_amount

    if debug:
        dirty_price = bond_values[maturity_step, 0]
        clean_price = dirty_price - accrued[maturity_step]
        print(m, _tree_times[m], accrued[m], dirty_price, clean_price, 0, 0)

    # Step back from maturity to expiry date but with no exercise allowed.
    for m in range(maturity_step - 1, expiry_step, -1):

        nm = min(m, j_max)
        flow = tree_flows[m] * face_amount

        for k in range(-nm, nm + 1):
            kn = k + j_max
            r = _rt[m, kn]
            df = np.exp(-r * _dt)

            pu = _pu[kn]
            pm = _pm[kn]
            pd = _pd[kn]

            if k == j_max:
                vu = bond_values[m + 1, kn]
                vm = bond_values[m + 1, kn - 1]
                vd = bond_values[m + 1, kn - 2]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v
            elif k == -j_max:
                vu = bond_values[m + 1, kn + 2]
                vm = bond_values[m + 1, kn + 1]
                vd = bond_values[m + 1, kn]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v
            else:
                vu = bond_values[m + 1, kn + 1]
                vm = bond_values[m + 1, kn]
                vd = bond_values[m + 1, kn - 1]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v

            bond_values[m, kn] += flow

        if debug:
            print(m, _tree_times[m], accrued[m], dirty_price, clean_price, 0, 0)

    # Now consider exercise of the option on and before the expiry date
    for m in range(expiry_step, -1, -1):
        nm = min(m, j_max)
        flow = tree_flows[m] * face_amount

        for k in range(-nm, nm + 1):

            kn = k + j_max
            r = _rt[m, kn]
            df = np.exp(-r * _dt)

            pu = _pu[kn]
            pm = _pm[kn]
            pd = _pd[kn]

            if k == j_max:
                vu = bond_values[m + 1, kn]
                vm = bond_values[m + 1, kn - 1]
                vd = bond_values[m + 1, kn - 2]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v
            elif k == -j_max:
                vu = bond_values[m + 1, kn + 2]
                vm = bond_values[m + 1, kn + 1]
                vd = bond_values[m + 1, kn]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v
            else:
                vu = bond_values[m + 1, kn + 1]
                vm = bond_values[m + 1, kn]
                vd = bond_values[m + 1, kn - 1]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, kn] = v

            bond_values[m, kn] += flow

            vcall = 0.0
            vput = 0.0

            if k == j_max:
                vu = call_option_values[m + 1, kn]
                vm = call_option_values[m + 1, kn - 1]
                vd = call_option_values[m + 1, kn - 2]
                vcall = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = call_option_values[m + 1, kn + 2]
                vm = call_option_values[m + 1, kn + 1]
                vd = call_option_values[m + 1, kn]
                vcall = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = call_option_values[m + 1, kn + 1]
                vm = call_option_values[m + 1, kn]
                vd = call_option_values[m + 1, kn - 1]
                vcall = (pu * vu + pm * vm + pd * vd) * df

            call_option_values[m, kn] = vcall

            if k == j_max:
                vu = put_option_values[m + 1, kn]
                vm = put_option_values[m + 1, kn - 1]
                vd = put_option_values[m + 1, kn - 2]
                vput = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = put_option_values[m + 1, kn + 2]
                vm = put_option_values[m + 1, kn + 1]
                vd = put_option_values[m + 1, kn]
                vput = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = put_option_values[m + 1, kn + 1]
                vm = put_option_values[m + 1, kn]
                vd = put_option_values[m + 1, kn - 1]
                vput = (pu * vu + pm * vm + pd * vd) * df

            put_option_values[m, kn] = vput

            dirty_price = bond_values[m, kn]
            clean_price = dirty_price - accrued[m]
            call_exercise = max(clean_price - strike_price, 0.0)
            put_exercise = max(strike_price - clean_price, 0.0)

            hold_call = call_option_values[m, kn]
            hold_put = put_option_values[m, kn]

            if m == expiry_step:

                call_option_values[m, kn] = max(call_exercise, hold_call)
                put_option_values[m, kn] = max(put_exercise, hold_put)

            elif exercise_type_int == 3 and m < expiry_step:  # AMERICAN

                call_option_values[m, kn] = max(call_exercise, hold_call)
                put_option_values[m, kn] = max(put_exercise, hold_put)

        if debug:
            print(
                m,
                _tree_times[m],
                accrued[m],
                dirty_price,
                clean_price,
                call_exercise,
                put_exercise,
            )

    return call_option_values[0, j_max], put_option_values[0, j_max]

########################################################################################

@njit(fastmath=True, cache=True)
def callable_puttable_bond_tree_fast(

    cpn_times,
    cpn_flows,
    call_times,
    call_prices,
    put_times,
    put_prices,
    face_amount,
    _sigma,
    _a,
    _qq,  # IS SIGMA USED ?
    _pu,
    _pm,
    _pd,
    _rt,
    _dt,
    _tree_times,
    _df_times,
    _df_values,
):
    """Value a bond with embedded put and call options that can be exercised
    at any time over the specified list of put and call dates.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time."""

    num_time_steps, num_nodes = _qq.shape
    dt = _dt
    j_max = ceil(0.1835 / (_a * _dt))
    t_mat = cpn_times[-1]
    maturity_step = int(t_mat / dt + 0.50)

    # Map coupons onto tree while preserving their present value

    tree_flows = np.zeros(num_time_steps)

    num_cpns = len(cpn_times)
    for i in range(0, num_cpns):
        t_cpn = cpn_times[i]
        n = int(t_cpn / _dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(t_cpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        tree_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree

    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent

    mapped_times = np.array([0.0])
    mapped_amounts = np.array([0.0])
    for n in range(1, len(_tree_times)):
        if tree_flows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, tree_flows[n])


    accrued = np.zeros(num_time_steps)
    for m in range(0, num_time_steps):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, mapped_times, mapped_amounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if tree_flows[m] > 0.0:
            accrued[m] = tree_flows[m] * face_amount

    # map call onto tree - must have no calls at high value

    tree_call_value = np.ones(num_time_steps) * face_amount * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(call_time / dt + 0.50)
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    tree_put_value = np.zeros(num_time_steps)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(put_time / dt + 0.50)
        tree_put_value[n] = put_prices[i]

    # Value the bond by backward induction starting at bond maturity

    call_put_bond_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    debug = False
    if debug:
        df = 1.0
        px = 0.0
        for i in range(0, maturity_step + 1):
            flow = tree_flows[i]
            t = _tree_times[i]
            df = _uinterpolate(t, _df_times, _df_values, interp)
            px += flow * df
        px += df

    # Now step back to today considering early exercise

    m = maturity_step
    nm = min(maturity_step, j_max)
    vcall = tree_call_value[m]
    vput = tree_put_value[m]
    vhold = (1.0 + tree_flows[m]) * face_amount
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm + 1):
        kn = k + j_max
        bond_values[m, kn] = (1.0 + tree_flows[m]) * face_amount
        call_put_bond_values[m, kn] = value

    for m in range(maturity_step - 1, -1, -1):
        nm = min(m, j_max)
        flow = tree_flows[m] * face_amount
        vcall = tree_call_value[m]
        vput = tree_put_value[m]

        for k in range(-nm, nm + 1):
            kn = k + j_max
            rt = _rt[m, kn]
            df = np.exp(-rt * dt)
            pu = _pu[kn]
            pm = _pm[kn]
            pd = _pd[kn]

            if k == j_max:
                vu = bond_values[m + 1, kn]
                vm = bond_values[m + 1, kn - 1]
                vd = bond_values[m + 1, kn - 2]
            elif k == -j_max:
                vu = bond_values[m + 1, kn + 2]
                vm = bond_values[m + 1, kn + 1]
                vd = bond_values[m + 1, kn]
            else:
                vu = bond_values[m + 1, kn + 1]
                vm = bond_values[m + 1, kn]
                vd = bond_values[m + 1, kn - 1]

            v = (pu * vu + pm * vm + pd * vd) * df
            bond_values[m, kn] = v
            bond_values[m, kn] += flow

            if k == j_max:
                vu = call_put_bond_values[m + 1, kn]
                vm = call_put_bond_values[m + 1, kn - 1]
                vd = call_put_bond_values[m + 1, kn - 2]
            elif k == -j_max:
                vu = call_put_bond_values[m + 1, kn + 2]
                vm = call_put_bond_values[m + 1, kn + 1]
                vd = call_put_bond_values[m + 1, kn]
            else:
                vu = call_put_bond_values[m + 1, kn + 1]
                vm = call_put_bond_values[m + 1, kn]
                vd = call_put_bond_values[m + 1, kn - 1]

            vhold = (pu * vu + pm * vm + pd * vd) * df
            # Need to make add on coupons paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            call_put_bond_values[m, kn] = value

    return {
        "bondwithoption": call_put_bond_values[0, j_max],
        "bondpure": bond_values[0, j_max],
    }

########################################################################################

@njit(fastmath=True, cache=True)
def build_tree_fast(a, sigma, tree_times, num_time_steps, discount_factors):

    """Calibrate the tree to a term structure of interest rates."""

    tree_maturity = tree_times[-1]
    dt = tree_maturity / (num_time_steps + 1)
    dx = sigma * np.sqrt(3.0 * dt)
    j_max = ceil(0.1835 / (a * dt))

    if j_max > 1000:
        raise FinError("j_max > 1000. Increase a or dt.")

    pu = np.zeros(shape=(2 * j_max + 1))
    pm = np.zeros(shape=(2 * j_max + 1))
    pd = np.zeros(shape=(2 * j_max + 1))

    # The short rate goes out one step extra to have the final short rate
    # This is the BK model so x = log(r)
    x = np.zeros(shape=(num_time_steps + 2, 2 * j_max + 1))
    rt = np.zeros(shape=(num_time_steps + 2, 2 * j_max + 1))

    # probabilities start at time 0 and go out to one step before T
    # Branching is simple trinomial out to time step m=1 after which
    # the top node and bottom node connect internally to two lower nodes
    # and two upper nodes respectively. The probabilities only depend on j

    for j in range(-j_max, j_max + 1):

        ajdt = a * j * dt
        jn = j + j_max

        if j == j_max:
            pu[jn] = 7.0 / 6.0 + 0.50 * (ajdt * ajdt - 3.0 * ajdt)
            pm[jn] = -1.0 / 3.0 - ajdt * ajdt + 2.0 * ajdt
            pd[jn] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt - ajdt)
        elif j == -j_max:
            pu[jn] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt + ajdt)
            pm[jn] = -1.0 / 3.0 - ajdt * ajdt - 2.0 * ajdt
            pd[jn] = 7.0 / 6.0 + 0.50 * (ajdt * ajdt + 3.0 * ajdt)
        else:
            pu[jn] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt - ajdt)
            pm[jn] = 2.0 / 3.0 - ajdt * ajdt
            pd[jn] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt + ajdt)

    # Arrow-Debreu array
    qq = np.zeros(shape=(num_time_steps + 2, 2 * j_max + 1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(num_time_steps + 1)

    # Time zero is trivial for the Arrow-Debreu price
    qq[0, j_max] = 1.0

    # Estimate short rate over first year
    r0 = -np.log(discount_factors[1]) / tree_times[1]

    # We initialise x0 with value of log of r0
    x0 = np.log(r0)

    # Big loop over time steps
    for m in range(0, num_time_steps + 1):

        nm = min(m, j_max)

        # Need to do drift adjustment which is non-linear and so requires
        # a root search algorithm to find value of x0.

        alpha[m] = search_root_deriv(
            x0, nm, qq[m], discount_factors[m + 1], dx, dt, j_max
        )

        x0 = alpha[m]

        for j in range(-nm, nm + 1):
            jn = j + j_max
            x[m, jn] = alpha[m] + j * dx
            rt[m, jn] = np.exp(x[m, jn])

        # Loop over all nodes at time m to calculate next values of qq
        for j in range(-nm, nm + 1):

            jn = j + j_max
            rdt = np.exp(x[m, jn]) * dt
            z = np.exp(-rdt)

            if j == j_max:
                qq[m + 1, jn] += qq[m, jn] * pu[jn] * z
                qq[m + 1, jn - 1] += qq[m, jn] * pm[jn] * z
                qq[m + 1, jn - 2] += qq[m, jn] * pd[jn] * z
            elif j == -j_max:
                qq[m + 1, jn] += qq[m, jn] * pd[jn] * z
                qq[m + 1, jn + 1] += qq[m, jn] * pm[jn] * z
                qq[m + 1, jn + 2] += qq[m, jn] * pu[jn] * z
            else:
                qq[m + 1, jn + 1] += qq[m, jn] * pu[jn] * z
                qq[m + 1, jn] += qq[m, jn] * pm[jn] * z
                qq[m + 1, jn - 1] += qq[m, jn] * pd[jn] * z

    return (qq, pu, pm, pd, rt, dt)

########################################################################################

class BKTree:

    ####################################################################################

    def __init__(self, sigma: float, a: float, num_time_steps: int = 100):

        """Constructs the Black Karasinski rate model. The speed of mean
        reversion a and volatility are passed in. The short rate process
        is given by d(log(r)) = (theta(t) - a*log(r)) * dt  + sigma * dW"""

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        if a < 1e-10:
            a = 1e-10

        self.a = a
        self.sigma = sigma

        if num_time_steps < 3:
            raise FinError("Drift fitting requires at least 3 time steps")

        self.num_time_steps = num_time_steps

        self.qq = None
        self.rt = None
        self.tree_times = None
        self.pu = None
        self.pm = None
        self.pd = None
        self.discount_curve = None
        self.df_times = None
        self.dfs = None
        self.dt = None

    ####################################################################################

    def build_tree(self, t_mat, df_times, df_values):

        if isinstance(df_times, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(df_values, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = InterpTypes.FLAT_FWD_RATES.value

        tree_maturity = t_mat * (self.num_time_steps + 1) / self.num_time_steps
        tree_times = np.linspace(0.0, tree_maturity, self.num_time_steps + 2)
        self.tree_times = tree_times

        df_tree = np.zeros(shape=(self.num_time_steps + 2))
        df_tree[0] = 1.0

        for i in range(1, self.num_time_steps + 2):
            t = tree_times[i]
            df_tree[i] = _uinterpolate(t, df_times, df_values, interp)

        self.df_times = df_times
        self.dfs = df_values

        self.qq, self.pu, self.pm, self.pd, self.rt, self.dt = build_tree_fast(
            self.a, self.sigma, tree_times, self.num_time_steps, df_tree
        )

        return

    ####################################################################################

    def bond_option(

        self,
        t_exp,
        strike_price,
        face_amount,
        cpn_times,
        cpn_flows,
        exercise_type,
    ):
        """Value a bond option that has European or American exercise using
        the Black-Karasinski model. The model uses a trinomial tree."""

        exercise_type_int = option_exercise_types_to_int(exercise_type)

        t_mat = cpn_times[-1]

        if t_exp > t_mat:
            raise FinError("Option expiry after bond matures.")

        if t_exp < 0.0:
            raise FinError("Option expiry time negative.")


        call_value, put_value = american_bond_option_tree_fast(
            t_exp,
            t_mat,
            strike_price,
            face_amount,
            cpn_times,
            cpn_flows,
            exercise_type_int,
            self.df_times,
            self.dfs,
            self.tree_times,
            self.qq,
            self.pu,
            self.pm,
            self.pd,
            self.rt,
            self.dt,
            self.a,
        )

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def bermudan_swaption(

        self,
        t_exp,
        t_mat,
        strike_price,
        face_amount,
        cpn_times,
        cpn_flows,
        exercise_type,
    ):
        """Swaption that can be exercised on specific dates over the exercise
        period. Due to non-analytical bond price we need to extend tree out to
        bond maturity and take into account cash flows through time."""

        exercise_type_int = option_exercise_types_to_int(exercise_type)

        t_mat = cpn_times[-1]

        if t_exp > t_mat:
            raise FinError("Option expiry after bond matures.")

        if t_exp < 0.0:
            raise FinError("Option expiry time negative.")


        pay_value, rec_value = bermudan_swaption_tree_fast(
            t_exp,
            t_mat,
            strike_price,
            face_amount,
            cpn_times,
            cpn_flows,
            exercise_type_int,
            self.df_times,
            self.dfs,
            self.tree_times,
            self.qq,
            self.pu,
            self.pm,
            self.pd,
            self.rt,
            self.dt,
            self.a,
        )

        return {"pay": pay_value, "rec": rec_value}

    ####################################################################################

    def callable_puttable_bond_tree(

        self,
        cpn_times,
        cpn_flows,
        call_times,
        call_prices,
        put_times,
        put_prices,
        face,
    ):
        """Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time."""

        call_times = np.array(call_times)
        put_times = np.array(put_times)

        call_prices = np.array(call_prices)
        put_prices = np.array(put_prices)

        v = callable_puttable_bond_tree_fast(
            cpn_times,
            cpn_flows,
            call_times,
            call_prices,
            put_times,
            put_prices,
            face,
            self.sigma,
            self.a,
            self.qq,
            self.pu,
            self.pm,
            self.pd,
            self.rt,
            self.dt,
            self.tree_times,
            self.df_times,
            self.dfs,
        )

        return {
            "bondwithoption": v["bondwithoption"],
            "bondpure": v["bondpure"],
        }

    ####################################################################################

    def __repr__(self):

        """Return string with class details."""

        s = "Black-Karasinski Model\n"
        s += label_to_string("Sigma", self.sigma)
        s += label_to_string("a", self.a)
        s += label_to_string("num_time_steps", self.num_time_steps)
        return s

