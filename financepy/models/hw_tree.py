# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from math import ceil

from enum import Enum
import numpy as np
from scipy import optimize
from numba import njit

from ..utils.error import FinError
from ..utils.math import normcdf, accrued_interpolator
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string
from ..utils.global_types import FinExerciseTypes
from ..utils.global_vars import G_SMALL

INTERP_TYPE_VALUE = InterpTypes.FLAT_FWD_RATES.value
SMALL = 1e-10

# TODO: Put Jamshidian code into Numba to get speed up

# dr = (theta(t) - r) dt + sigma * dW

########################################################################################


class FinHWEuropeanCalcType(Enum):

    JAMSHIDIAN = 1
    EXPIRY_ONLY = 2
    EXPIRY_TREE = 3


########################################################################################


def option_exercise_types_to_int(option_exercise_type):

    if option_exercise_type == FinExerciseTypes.EUROPEAN:
        return 1
    if option_exercise_type == FinExerciseTypes.BERMUDAN:
        return 2
    if option_exercise_type == FinExerciseTypes.AMERICAN:
        return 3

    raise FinError("Unknown option exercise type.")


########################################################################################


@njit(fastmath=True, cache=True)
def p_fast(t, t_mat, r_t, delta, pt, ptd, p_mat, _sigma, _a):
    """Forward discount factor as seen at some time t which may be in the
    future for payment at time T where r_t is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor."""

    if abs(_a) < SMALL:
        _a = SMALL

    bt_t = (1.0 - np.exp(-_a * (t_mat - t))) / _a
    bt_delta = (1.0 - np.exp(-_a * delta)) / _a

    term1 = np.log(p_mat / pt) - (bt_t / bt_delta) * np.log(ptd / pt)

    term2 = (
        (_sigma**2) * (1 - np.exp(-2 * _a * t)) * bt_t * (bt_t - bt_delta) / (4 * _a)
    )

    log_ahat = term1 - term2
    bhatt_t = (bt_t / bt_delta) * delta
    p = np.exp(log_ahat - bhatt_t * r_t)
    return p


########################################################################################


@njit(fastmath=True, cache=True)
def build_tree_fast(a, sigma, tree_times, num_time_steps, discount_factors):
    """Fast tree construction using Numba."""
    tree_maturity = tree_times[-1]
    dt = tree_maturity / (num_time_steps + 1)
    d_r = sigma * np.sqrt(3.0 * dt)
    j_max = ceil(0.1835 / (a * dt))
    n = j_max

    pu = np.zeros(shape=2 * j_max + 1)
    pm = np.zeros(shape=2 * j_max + 1)
    pd = np.zeros(shape=2 * j_max + 1)

    # The short rate goes out one step extra to have the final short rate
    r_t = np.zeros(shape=(num_time_steps + 2, 2 * j_max + 1))

    # probabilities star_t at time 0 and go out to one step before T
    # Branching is simple trinomial out to time step m=1 after which
    # the top node and bottom node connect internally to two lower nodes
    # and two upper nodes respectively. The probabilities only depend on j

    for j in range(-j_max, j_max + 1):
        ajdt = a * j * dt
        j_n = j + n
        if j == j_max:
            pu[j_n] = 7.0 / 6.0 + 0.50 * (ajdt * ajdt - 3.0 * ajdt)
            pm[j_n] = -1.0 / 3.0 - ajdt * ajdt + 2.0 * ajdt
            pd[j_n] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt - ajdt)
        elif j == -j_max:
            pu[j_n] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt + ajdt)
            pm[j_n] = -1.0 / 3.0 - ajdt * ajdt - 2.0 * ajdt
            pd[j_n] = 7.0 / 6.0 + 0.50 * (ajdt * ajdt + 3.0 * ajdt)
        else:
            pu[j_n] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt - ajdt)
            pm[j_n] = 2.0 / 3.0 - ajdt * ajdt
            pd[j_n] = 1.0 / 6.0 + 0.50 * (ajdt * ajdt + ajdt)

    # Arrow-Debreu array
    qq = np.zeros(shape=(num_time_steps + 2, 2 * n + 1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(num_time_steps + 1)

    # Time zero is trivial for the Arrow-Debreu price
    qq[0, n] = 1.0

    # Big loop over time steps
    for m in range(0, num_time_steps + 1):

        nm = min(m, j_max)
        sum_qz = 0.0
        for j in range(-nm, nm + 1):
            rdt = j * d_r * dt
            sum_qz += qq[m, j + n] * np.exp(-rdt)
        alpha[m] = np.log(sum_qz / discount_factors[m + 1]) / dt

        for j in range(-nm, nm + 1):
            j_n = j + n
            r_t[m, j_n] = alpha[m] + j * d_r

        # Loop over all nodes at time m to calculate next values of qq
        for j in range(-nm, nm + 1):
            j_n = j + n
            rdt = r_t[m, j_n] * dt
            z = np.exp(-rdt)

            if j == j_max:
                qq[m + 1, j_n] += qq[m, j_n] * pu[j_n] * z
                qq[m + 1, j_n - 1] += qq[m, j_n] * pm[j_n] * z
                qq[m + 1, j_n - 2] += qq[m, j_n] * pd[j_n] * z
            elif j == -j_max:
                qq[m + 1, j_n] += qq[m, j_n] * pd[j_n] * z
                qq[m + 1, j_n + 1] += qq[m, j_n] * pm[j_n] * z
                qq[m + 1, j_n + 2] += qq[m, j_n] * pu[j_n] * z
            else:
                qq[m + 1, j_n + 1] += qq[m, j_n] * pu[j_n] * z
                qq[m + 1, j_n] += qq[m, j_n] * pm[j_n] * z
                qq[m + 1, j_n - 1] += qq[m, j_n] * pd[j_n] * z

    return (qq, pu, pm, pd, r_t, dt)


########################################################################################


@njit(fastmath=True, cache=True)
def american_bond_option_tree_fast(
    t_exp,
    strike_price,
    face_amount,
    cpn_times,
    cpn_amounts,
    exercise_type_int,
    _sigma,
    _a,
    _qq,
    _pu,
    _pm,
    _pd,
    _r_t,
    _dt,
    _tree_times,
    _df_times,
    _df_values,
):
    """Value an option on a bond with cpns that can have European or
    American exercise. Some minor issues to do with handling cpns on
    the option expiry date need to be solved."""

    debug = False
    if debug:
        print("Entering AmerBondOption")
        print("Coupon Times", cpn_times)
        print("Coupon Amounts", cpn_amounts)

    num_time_steps, num_nodes = _qq.shape
    dt = _dt
    j_max = ceil(0.1835 / (_a * dt))
    expiry_step = int(t_exp / dt + 0.50)

    # Want to add cpns before expiry to the grid so that we can value
    # their impact on the decision to exercise the option early
    tree_flows = np.zeros(num_time_steps)
    num_cpns = len(cpn_times)

    # Flows that fall on the expiry date included. The tree only goes out to
    # the expiry date so cpns after this date do not go onto the tree.
    for i in range(0, num_cpns):
        t_cpn = cpn_times[i]
        if t_cpn <= t_exp:
            n = int(t_cpn / dt + 0.50)
            ttree = _tree_times[n]
            df_flow = _uinterpolate(t_cpn, _df_times, _df_values, INTERP_TYPE_VALUE)
            df_tree = _uinterpolate(ttree, _df_times, _df_values, INTERP_TYPE_VALUE)
            tree_flows[n] += cpn_amounts[i] * 1.0 * df_flow / df_tree

    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent

    # I star_t the tree with the previous cpn time and amount
    # (does not matter)
    mapped_times = np.zeros(0)  # CHANGE
    mapped_amounts = np.zeros(0)  # CHANGE
    for n in range(0, len(_tree_times)):
        if tree_flows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, tree_flows[n])

    # Need future cash flows which are exact time and size for accrued at t_exp
    for n in range(0, num_cpns):
        if cpn_times[n] > t_exp:
            mapped_times = np.append(mapped_times, cpn_times[n])
            mapped_amounts = np.append(mapped_amounts, cpn_amounts[n])

    if debug:
        print("MAPPED TIMES", mapped_times)
        print("MAPPED AMOUNTS", mapped_amounts)

    accrued = np.zeros(num_time_steps)
    for m in range(0, expiry_step + 1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, cpn_times, cpn_amounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if tree_flows[m] > 0.0:
            accrued[m] = tree_flows[m] * face_amount

    if debug:
        for i in range(0, expiry_step + 1):
            print(i, tree_flows[i], accrued[i])

    call_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    put_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    pt_exp = _uinterpolate(t_exp, _df_times, _df_values, INTERP_TYPE_VALUE)
    ptdelta = _uinterpolate(t_exp + dt, _df_times, _df_values, INTERP_TYPE_VALUE)

    cpn = 0.0
    zcb = 0.0

    # As the HW model has a closed form solution for the bond price, I use
    # this fact to calculate the bond price at expiry on the tree nodes

    nm = min(expiry_step, j_max)
    for k in range(-nm, nm + 1):
        k_n = k + j_max
        r_t = _r_t[expiry_step, k_n]
        bond_price = 0.0
        for i in range(0, num_cpns):
            tflow = cpn_times[i]
            if tflow >= t_exp:
                ptflow = _uinterpolate(tflow, _df_times, _df_values, INTERP_TYPE_VALUE)
                zcb = p_fast(t_exp, tflow, r_t, dt, pt_exp, ptdelta, ptflow, _sigma, _a)
                cpn = cpn_amounts[i]
                bond_price += cpn * face_amount * zcb

        bond_price += zcb * face_amount

        # The flow on this date has been added
        bond_values[expiry_step, k_n] = bond_price

    # Now consider exercise of the option on the expiry date
    nm = min(expiry_step, j_max)
    for k in range(-nm, nm + 1):
        k_n = k + j_max
        dirty_price = bond_values[expiry_step, k_n]
        clean_price = dirty_price - accrued[expiry_step]
        call_exercise = max(clean_price - strike_price, 0.0)
        put_exercise = max(strike_price - clean_price, 0.0)
        call_option_values[expiry_step, k_n] = call_exercise
        put_option_values[expiry_step, k_n] = put_exercise

    m = expiry_step

    if debug:
        print("-----------------------------------------")
        print(
            "EXP",
            _tree_times[m],
            accrued[m],
            dirty_price,
            clean_price,
            call_exercise,
            put_exercise,
        )

    #        print(k_n, bond_values[expiry_step, k_n], "CLEAN", clean_price)
    #        print("EXPIRY DATE", k_n, clean_price, accrued[expiry_step], strike_price)

    # Now step back to today considering exercise at expiry and before
    for m in range(expiry_step - 1, -1, -1):
        nm = min(m, j_max)
        flow = tree_flows[m] * face_amount

        for k in range(-nm, nm + 1):
            k_n = k + j_max
            r = _r_t[m, k_n]
            df = np.exp(-r * dt)

            pu = _pu[k_n]
            pm = _pm[k_n]
            pd = _pd[k_n]

            if k == j_max:
                vu = bond_values[m + 1, k_n]
                vm = bond_values[m + 1, k_n - 1]
                vd = bond_values[m + 1, k_n - 2]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, k_n] = v
            elif k == -j_max:
                vu = bond_values[m + 1, k_n + 2]
                vm = bond_values[m + 1, k_n + 1]
                vd = bond_values[m + 1, k_n]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, k_n] = v
            else:
                vu = bond_values[m + 1, k_n + 1]
                vm = bond_values[m + 1, k_n]
                vd = bond_values[m + 1, k_n - 1]
                v = (pu * vu + pm * vm + pd * vd) * df
                bond_values[m, k_n] = v

            bond_values[m, k_n] += flow

            vcall = 0.0
            vput = 0.0

            if k == j_max:
                vu = call_option_values[m + 1, k_n]
                vm = call_option_values[m + 1, k_n - 1]
                vd = call_option_values[m + 1, k_n - 2]
                vcall = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = call_option_values[m + 1, k_n + 2]
                vm = call_option_values[m + 1, k_n + 1]
                vd = call_option_values[m + 1, k_n]
                vcall = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = call_option_values[m + 1, k_n + 1]
                vm = call_option_values[m + 1, k_n]
                vd = call_option_values[m + 1, k_n - 1]
                vcall = (pu * vu + pm * vm + pd * vd) * df

            call_option_values[m, k_n] = vcall

            if k == j_max:
                vu = put_option_values[m + 1, k_n]
                vm = put_option_values[m + 1, k_n - 1]
                vd = put_option_values[m + 1, k_n - 2]
                vput = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = put_option_values[m + 1, k_n + 2]
                vm = put_option_values[m + 1, k_n + 1]
                vd = put_option_values[m + 1, k_n]
                vput = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = put_option_values[m + 1, k_n + 1]
                vm = put_option_values[m + 1, k_n]
                vd = put_option_values[m + 1, k_n - 1]
                vput = (pu * vu + pm * vm + pd * vd) * df

            put_option_values[m, k_n] = vput

            dirty_price = bond_values[m, k_n]
            clean_price = dirty_price - accrued[m]
            call_exercise = max(clean_price - strike_price, 0.0)
            put_exercise = max(strike_price - clean_price, 0.0)

            hold_call = call_option_values[m, k_n]
            hold_put = put_option_values[m, k_n]

            if m == expiry_step:

                call_option_values[m, k_n] = max(call_exercise, hold_call)
                put_option_values[m, k_n] = max(put_exercise, hold_put)

            elif exercise_type_int == 3 and m < expiry_step:  # AMERICAN

                call_option_values[m, k_n] = max(call_exercise, hold_call)
                put_option_values[m, k_n] = max(put_exercise, hold_put)

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
    _r_t,
    _dt,
    _a,
):
    """Option to enter into a swap that can be exercised on cpn payment
    dates after the star_t of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time."""

    num_time_steps, num_nodes = _qq.shape
    j_max = ceil(0.1835 / (_a * _dt))
    expiry_step = int(t_exp / _dt + 0.50)
    maturity_step = int(t_mat / _dt + 0.50)

    fixed_leg_flows = np.zeros(num_time_steps)
    float_leg_values = np.zeros(num_time_steps)
    num_cpns = len(cpn_times)

    # Tree flows go all the way out to the bond maturity date
    for i in range(0, num_cpns):
        t_cpn = cpn_times[i]
        n = int(round(t_cpn / _dt, 0))
        ttree = _tree_times[n]
        df_flow = _uinterpolate(t_cpn, _df_times, _df_values, INTERP_TYPE_VALUE)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, INTERP_TYPE_VALUE)
        fixed_leg_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree
        float_leg_values[n] = strike_price * df_flow / df_tree

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

    # The value of the swap at each time and node. Principal is exchanged.
    fixed_leg_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a payer swap
    pay_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a receiver swap
    rec_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Star_t with the value of the bond at maturity
    for k in range(0, num_nodes):
        flow = 1.0 + fixed_leg_flows[maturity_step]
        fixed_leg_values[maturity_step, k] = flow * face_amount

    n = j_max

    # Now step back to today considering early exercise
    for m in range(maturity_step - 1, -1, -1):
        nm = min(m, j_max)
        flow = fixed_leg_flows[m] * face_amount

        for k in range(-nm, nm + 1):
            k_n = k + n
            r_t = _r_t[m, k_n]
            df = np.exp(-r_t * _dt)
            pu = _pu[k_n]
            pm = _pm[k_n]
            pd = _pd[k_n]

            if k == j_max:
                vu = fixed_leg_values[m + 1, k_n]
                vm = fixed_leg_values[m + 1, k_n - 1]
                vd = fixed_leg_values[m + 1, k_n - 2]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, k_n] = v
            elif k == -j_max:
                vu = fixed_leg_values[m + 1, k_n + 2]
                vm = fixed_leg_values[m + 1, k_n + 1]
                vd = fixed_leg_values[m + 1, k_n]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, k_n] = v
            else:
                vu = fixed_leg_values[m + 1, k_n + 1]
                vm = fixed_leg_values[m + 1, k_n]
                vd = fixed_leg_values[m + 1, k_n - 1]
                v = (pu * vu + pm * vm + pd * vd) * df
                fixed_leg_values[m, k_n] = v

            fixed_leg_values[m, k_n] += flow
            vpay = 0.0
            vrec = 0.0

            if k == j_max:
                vu = pay_values[m + 1, k_n]
                vm = pay_values[m + 1, k_n - 1]
                vd = pay_values[m + 1, k_n - 2]
                vpay = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = pay_values[m + 1, k_n + 2]
                vm = pay_values[m + 1, k_n + 1]
                vd = pay_values[m + 1, k_n]
                vpay = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = pay_values[m + 1, k_n + 1]
                vm = pay_values[m + 1, k_n]
                vd = pay_values[m + 1, k_n - 1]
                vpay = (pu * vu + pm * vm + pd * vd) * df

            pay_values[m, k_n] = vpay

            if k == j_max:
                vu = rec_values[m + 1, k_n]
                vm = rec_values[m + 1, k_n - 1]
                vd = rec_values[m + 1, k_n - 2]
                vrec = (pu * vu + pm * vm + pd * vd) * df
            elif k == -j_max:
                vu = rec_values[m + 1, k_n + 2]
                vm = rec_values[m + 1, k_n + 1]
                vd = rec_values[m + 1, k_n]
                vrec = (pu * vu + pm * vm + pd * vd) * df
            else:
                vu = rec_values[m + 1, k_n + 1]
                vm = rec_values[m + 1, k_n]
                vd = rec_values[m + 1, k_n - 1]
                vrec = (pu * vu + pm * vm + pd * vd) * df

            rec_values[m, k_n] = vrec

            hold_pay = pay_values[m, k_n]
            hold_rec = rec_values[m, k_n]

            # The floating value is clean and so must be the fixed value
            fixed_leg_value = fixed_leg_values[m, k_n] - accrued[m]
            float_leg_value = float_leg_values[m]

            pay_exercise = max(float_leg_value - fixed_leg_value, 0.0)
            rec_exercise = max(fixed_leg_value - float_leg_value, 0.0)

            if m == expiry_step:

                pay_values[m, k_n] = max(pay_exercise, hold_pay)
                rec_values[m, k_n] = max(rec_exercise, hold_rec)

            elif exercise_type_int == 2 and flow > G_SMALL and m >= expiry_step:

                pay_values[m, k_n] = max(pay_exercise, hold_pay)
                rec_values[m, k_n] = max(rec_exercise, hold_rec)

            elif exercise_type_int == 3 and m >= expiry_step:

                pay_values[m, k_n] = max(pay_exercise, hold_pay)
                rec_values[m, k_n] = max(rec_exercise, hold_rec)

                # Need to define floating value on all grid dates

                raise FinError("American optionality not tested.")

    return pay_values[0, j_max], rec_values[0, j_max]


# TODO: CHECK ACCRUED AND COUPONS TO SEE IF IT WORKS FOR LOW TREE STEPS

########################################################################################


@njit(fastmath=True, cache=True)
def callable_puttable_bond_tree_fast(
    cpn_times,
    cpn_flows,
    call_times,
    call_prices,
    put_times,
    put_prices,
    face,
    _sigma,
    _a,
    _qq,  # IS SIGMA USED ?
    _pu,
    _pm,
    _pd,
    _r_t,
    _dt,
    _tree_times,
    _df_times,
    _df_values,
):
    """Value an option on a bond with cpns that can have European or
    American exercise. Some minor issues to do with handling cpns on
    the option expiry date need to be solved."""

    #    print("Coupon Times:", cpn_times)
    #    print("Coupon Flows:", cpn_flows)

    #    print("DF Times:", _df_times)
    #    print("DF Values:", _df_values)

    if np.any(cpn_times < 0.0):
        raise FinError("No cpn times can be before the value date.")

    num_time_steps, num_nodes = _qq.shape
    dt = _dt
    j_max = ceil(0.1835 / (_a * dt))
    t_mat = cpn_times[-1]
    maturity_step = int(t_mat / dt + 0.50)

    # Map cpns onto tree while preserving their present value

    tree_flows = np.zeros(num_time_steps)

    num_cpns = len(cpn_times)
    for i in range(0, num_cpns):
        t_cpn = cpn_times[i]
        n = int(round(t_cpn / dt, 0))
        ttree = _tree_times[n]
        df_flow = _uinterpolate(t_cpn, _df_times, _df_values, INTERP_TYPE_VALUE)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, INTERP_TYPE_VALUE)
        tree_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree

    #    print("Tree flows:", tree_flows)

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
        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if tree_flows[m] > 0.0:
            accrued[m] = tree_flows[m] * face

    # map call onto tree - must have no calls at high value

    tree_call_value = np.ones(num_time_steps) * face * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(round(call_time / dt, 0))
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    tree_put_value = np.zeros(num_time_steps)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(round(put_time / dt, 0))
        tree_put_value[n] = put_prices[i]

    # Value the bond by backward induction star_ting at bond maturity

    call_put_bond_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    debug = True
    if debug:
        df = 1.0
        px = 0.0
        for i in range(0, maturity_step + 1):
            flow = tree_flows[i]
            t = _tree_times[i]
            df = _uinterpolate(t, _df_times, _df_values, INTERP_TYPE_VALUE)

            if flow > G_SMALL:
                pv = flow * df
                px += pv

        px += df

    # Now step back to today considering early exercise

    m = maturity_step
    nm = min(maturity_step, j_max)
    vcall = tree_call_value[m]
    vput = tree_put_value[m]
    vhold = (1.0 + tree_flows[m]) * face
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm + 1):
        k_n = k + j_max
        bond_values[m, k_n] = (1.0 + tree_flows[m]) * face
        call_put_bond_values[m, k_n] = value

    # Now step back to today considering early put and call
    for m in range(maturity_step - 1, -1, -1):
        nm = min(m, j_max)
        flow = tree_flows[m] * face
        vcall = tree_call_value[m]
        vput = tree_put_value[m]

        for k in range(-nm, nm + 1):
            k_n = k + j_max
            r_t = _r_t[m, k_n]
            df = np.exp(-r_t * dt)
            pu = _pu[k_n]
            pm = _pm[k_n]
            pd = _pd[k_n]

            if k == j_max:
                vu = bond_values[m + 1, k_n]
                vm = bond_values[m + 1, k_n - 1]
                vd = bond_values[m + 1, k_n - 2]
            elif k == -j_max:
                vu = bond_values[m + 1, k_n + 2]
                vm = bond_values[m + 1, k_n + 1]
                vd = bond_values[m + 1, k_n]
            else:
                vu = bond_values[m + 1, k_n + 1]
                vm = bond_values[m + 1, k_n]
                vd = bond_values[m + 1, k_n - 1]

            v = (pu * vu + pm * vm + pd * vd) * df
            bond_values[m, k_n] = v
            bond_values[m, k_n] += flow

            if k == j_max:
                vu = call_put_bond_values[m + 1, k_n]
                vm = call_put_bond_values[m + 1, k_n - 1]
                vd = call_put_bond_values[m + 1, k_n - 2]
            elif k == -j_max:
                vu = call_put_bond_values[m + 1, k_n + 2]
                vm = call_put_bond_values[m + 1, k_n + 1]
                vd = call_put_bond_values[m + 1, k_n]
            else:
                vu = call_put_bond_values[m + 1, k_n + 1]
                vm = call_put_bond_values[m + 1, k_n]
                vd = call_put_bond_values[m + 1, k_n - 1]

            vhold = (pu * vu + pm * vm + pd * vd) * df
            # Need to make add on cpns paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            call_put_bond_values[m, k_n] = value

    return {
        "bondwithoption": call_put_bond_values[0, j_max],
        "bondpure": bond_values[0, j_max],
    }


########################################################################################


def fwd_dirty_bond_price(r_t, *args):
    """Price a cpn bearing bond on the option expiry date and return
    the difference from a strike price. This is used in a root search to
    find the future expiry time short rate that makes the bond price equal
    to the option strike price. It is a key step in the Jamshidian bond
    decomposition approach. The strike is a clean price."""

    self = args[0]
    t_exp = args[1]
    cpn_times = args[2]
    cpn_amounts = args[3]
    df_times = args[4]
    df_values = args[5]
    strike_price = args[6]
    face = args[7]

    dt = 0.001
    tdelta = t_exp + dt
    pt_exp = _uinterpolate(t_exp, df_times, df_values, INTERP_TYPE_VALUE)
    ptdelta = _uinterpolate(tdelta, df_times, df_values, INTERP_TYPE_VALUE)

    #    print("TEXP", t_exp, pt_exp)

    num_flows = len(cpn_times)
    pv = 0.0

    for i in range(1, num_flows):

        t_cpn = cpn_times[i]
        cpn = cpn_amounts[i]

        if t_cpn > t_exp:
            pt_cpn = _uinterpolate(t_cpn, df_times, df_values, INTERP_TYPE_VALUE)
            zcb = p_fast(
                t_exp,
                t_cpn,
                r_t,
                dt,
                pt_exp,
                ptdelta,
                pt_cpn,
                self.sigma,
                self.a,
            )
            pv = pv + zcb * cpn
    #            print("TCPN", t_cpn, "ZCB", zcb, "CPN", cpn, "PV", pv)

    if t_cpn >= t_exp:
        pv = pv + zcb

    #    print("TCPN", t_cpn, "ZCB", zcb, "PRI", 1.0, "PV", pv)

    accd = accrued_interpolator(t_exp, cpn_times, cpn_amounts)
    #    print("Accrued:", accd)
    #    print("t_exp:", t_exp)
    #    print("cpn_times:", cpn_times)
    #    print("cpn_amounts:", cpn_amounts)

    pv_clean = pv - accd
    obj = face * pv_clean - strike_price

    #    print("FWD PRICE", r_t, pv, accd, strike_price, obj)
    return obj


########################################################################################


class HWTree:

    ####################################################################################

    def __init__(
        self,
        sigma,
        a,
        num_time_steps=100,
        european_calc_type=FinHWEuropeanCalcType.EXPIRY_TREE,
    ):
        """Constructs the Hull-White rate model. The speed of mean reversion
        a and volatility are passed in. The short rate process is given by
        dr = (theta(t) - ar) * dt  + sigma * dW. The model will switch to use
        Jamshidian's approach where possible unless the useJamshidian flag is
        set to false in which case it uses the trinomial Tree."""

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        self.sigma = sigma
        self.a = a
        self.num_time_steps = num_time_steps
        self.european_calc_type = european_calc_type

        self.qq = None
        self.r = None
        self.tree_times = None
        self.pu = None
        self.pm = None
        self.pd = None
        self.discount_curve = None
        self.tree_built = False
        self.df_times = None
        self.dfs = None
        self.r_t = None
        self.dt = None

    ####################################################################################

    def option_on_zcb(self, t_exp, t_mat, strike, face_amount, df_times, df_values):
        """Price an option on a zero cpn bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date."""

        if t_exp > t_mat:
            raise FinError("Option expiry after bond matures.")

        if t_exp < 0.0:
            raise FinError("Option expiry time negative.")

        pt_exp = _uinterpolate(t_exp, df_times, df_values, INTERP_TYPE_VALUE)
        pt_mat = _uinterpolate(t_mat, df_times, df_values, INTERP_TYPE_VALUE)

        sigma = self.sigma
        a = self.a

        if abs(a) < SMALL:
            a = SMALL

        sigmap = (sigma / a) * (1.0 - np.exp(-a * (t_mat - t_exp)))
        sigmap = sigmap * np.sqrt((1.0 - np.exp(-2.0 * a * t_exp)) / 2.0 / a)

        if abs(sigmap) < SMALL:
            sigmap = SMALL

        h = np.log((face_amount * pt_mat) / (strike * pt_exp)) / sigmap + sigmap / 2.0
        call_value = face_amount * pt_mat * normcdf(h) - strike * pt_exp * normcdf(
            h - sigmap
        )
        put_value = strike * pt_exp * normcdf(
            -h + sigmap
        ) - face_amount * pt_mat * normcdf(-h)

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def european_bond_option_jamshidian(
        self, t_exp, strike_price, face, cpn_times, cpn_amounts, df_times, df_values
    ):
        """Valuation of a European bond option using the Jamshidian
        deconstruction of the bond into a strip of zero cpn bonds with the
        short rate that would make the bond option be at the money forward."""

        #        print(df_times)
        #        print(df_values)

        num_cpns = len(cpn_times)

        argtuple = (
            self,
            t_exp,
            cpn_times,
            cpn_amounts,
            df_times,
            df_values,
            strike_price,
            face,
        )

        # Can I improve on this initial guess ?
        x0 = 0.05

        rstar = optimize.newton(
            fwd_dirty_bond_price,
            x0=x0,
            fprime=None,
            args=argtuple,
            tol=1e-10,
            maxiter=50,
            fprime2=None,
        )

        # Now we price a series of zero cpn bonds using this short rate
        dt = 1e-6

        pt_exp = _uinterpolate(t_exp, df_times, df_values, INTERP_TYPE_VALUE)
        ptdelta = _uinterpolate(t_exp + dt, df_times, df_values, INTERP_TYPE_VALUE)

        call_value = 0.0
        put_value = 0.0

        # Adjust strike to handle
        for i in range(0, num_cpns):

            t_cpn = cpn_times[i]
            cpn = cpn_amounts[i]

            if t_cpn >= t_exp:  # cpns on the expiry date are included

                pt_cpn = _uinterpolate(t_cpn, df_times, df_values, INTERP_TYPE_VALUE)

                strike = p_fast(
                    t_exp,
                    t_cpn,
                    rstar,
                    dt,
                    pt_exp,
                    ptdelta,
                    pt_cpn,
                    self.sigma,
                    self.a,
                )

                v = self.option_on_zcb(t_exp, t_cpn, strike, 1.0, df_times, df_values)

                call = v["call"]
                put = v["put"]

                call_value += call * cpn * face
                put_value += put * cpn * face

        call_value += call * face
        put_value += put * face

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def european_bond_option_expiry_only(
        self, t_exp, strike_price, face_amount, cpn_times, cpn_amounts
    ):
        """Price a European option on a cpn-paying bond using a tree to
        generate short rates at the expiry date and then to use the analytical
        solution of zero cpn bond prices in the HW model to calculate the
        corresponding bond price. User provides bond object and option details.
        """

        dt = self.dt
        tdelta = t_exp + dt

        pt_exp = _uinterpolate(t_exp, self.df_times, self.dfs, INTERP_TYPE_VALUE)
        ptdelta = _uinterpolate(tdelta, self.df_times, self.dfs, INTERP_TYPE_VALUE)

        _, num_nodes = self.qq.shape
        expiry_step = int(t_exp / dt + 0.50)

        call_value = 0.0
        put_value = 0.0
        num_cpns = len(cpn_times)

        for k in range(0, num_nodes):

            q = self.qq[expiry_step, k]
            r_t = self.r_t[expiry_step, k]

            pv = 0.0

            for i in range(0, num_cpns):

                t_cpn = cpn_times[i]
                cpn = cpn_amounts[i]

                if t_cpn >= t_exp:

                    pt_cpn = _uinterpolate(
                        t_cpn, self.df_times, self.dfs, INTERP_TYPE_VALUE
                    )

                    zcb = p_fast(
                        t_exp,
                        t_cpn,
                        r_t,
                        dt,
                        pt_exp,
                        ptdelta,
                        pt_cpn,
                        self.sigma,
                        self.a,
                    )

                    pv += cpn * zcb

            pv += zcb

            #            print(t_exp)
            #            print(cpn_times)
            #            print(cpn_amounts)

            accrued = accrued_interpolator(t_exp, cpn_times, cpn_amounts)

            pv = pv - accrued
            #            print(accrued)

            put_payoff = max(strike_price - pv * face_amount, 0.0)
            call_payoff = max(pv * face_amount - strike_price, 0.0)

            put_value += q * put_payoff
            call_value += q * call_payoff

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def option_on_zero_cpn_bond_tree(self, t_exp, t_mat, strike_price, face_amount):
        """Price an option on a zero cpn bond using a HW trinomial
        tree. The discount curve was already supplied to the tree build."""

        if t_exp > t_mat:
            raise FinError("Option expiry after bond matures.")

        if t_exp < 0.0:
            raise FinError("Option expiry time negative.")

        if self.tree_times is None:
            raise FinError("Tree has not been constructed.")

        if self.tree_times[-1] < t_exp:
            raise FinError("Tree expiry must be >= option expiry date.")

        dt = self.dt
        tdelta = t_exp + dt

        pt_exp = _uinterpolate(t_exp, self.df_times, self.dfs, INTERP_TYPE_VALUE)
        ptdelta = _uinterpolate(tdelta, self.df_times, self.dfs, INTERP_TYPE_VALUE)
        pt_mat = _uinterpolate(t_mat, self.df_times, self.dfs, INTERP_TYPE_VALUE)

        _, num_nodes = self.qq.shape
        expiry_step = int(t_exp / dt + 0.50)

        call_value = 0.0
        put_value = 0.0

        for k in range(0, num_nodes):

            q = self.qq[expiry_step, k]
            r_t = self.r_t[expiry_step, k]

            zcb = p_fast(
                t_exp,
                t_mat,
                r_t,
                dt,
                pt_exp,
                ptdelta,
                pt_mat,
                self.sigma,
                self.a,
            )

            put_payoff = max(strike_price - zcb * face_amount, 0.0)
            call_payoff = max(zcb * face_amount - strike_price, 0.0)
            put_value += q * put_payoff
            call_value += q * call_payoff

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def bermudan_swaption(
        self, t_exp, t_mat, strike, face, cpn_times, cpn_flows, exercise_type
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
            strike,
            face,
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
            self.r_t,
            self.dt,
            self.a,
        )

        return {"pay": pay_value, "rec": rec_value}

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
        """Value a bond option that can have European or American exercise.
        This is done using a trinomial tree that we extend out to bond
        maturity. For European bond options, Jamshidian's model is
        faster and is used instead i.e. not this function."""

        exercise_type_int = option_exercise_types_to_int(exercise_type)

        if exercise_type_int == 1:

            if self.european_calc_type == FinHWEuropeanCalcType.JAMSHIDIAN:

                v = self.european_bond_option_jamshidian(
                    t_exp,
                    strike_price,
                    face_amount,
                    cpn_times,
                    cpn_flows,
                    self.df_times,
                    self.dfs,
                )

                call_value = v["call"]
                put_value = v["put"]

            elif self.european_calc_type == FinHWEuropeanCalcType.EXPIRY_ONLY:

                v = self.european_bond_option_expiry_only(
                    t_exp, strike_price, face_amount, cpn_times, cpn_flows
                )

                call_value = v["call"]
                put_value = v["put"]

            elif self.european_calc_type == FinHWEuropeanCalcType.EXPIRY_TREE:

                call_value, put_value = american_bond_option_tree_fast(
                    t_exp,
                    strike_price,
                    face_amount,
                    cpn_times,
                    cpn_flows,
                    exercise_type_int,
                    self.sigma,
                    self.a,
                    self.qq,
                    self.pu,
                    self.pm,
                    self.pd,
                    self.r_t,
                    self.dt,
                    self.tree_times,
                    self.df_times,
                    self.dfs,
                )

            else:
                raise FinError("Unknown HW model implementation choice.")

        else:

            call_value, put_value = american_bond_option_tree_fast(
                t_exp,
                strike_price,
                face_amount,
                cpn_times,
                cpn_flows,
                exercise_type_int,
                self.sigma,
                self.a,
                self.qq,
                self.pu,
                self.pm,
                self.pd,
                self.r_t,
                self.dt,
                self.tree_times,
                self.df_times,
                self.dfs,
            )

        return {"call": call_value, "put": put_value}

    ####################################################################################

    def callable_puttable_bond_tree(
        self,
        cpn_times,
        cpn_flows,
        call_times,
        call_prices,
        put_times,
        put_prices,
        face_amount,
    ):
        """Value an option on a bond with cpns that can have European or
        American exercise. Some minor issues to do with handling cpns on
        the option expiry date need to be solved. Also this function should be
        moved out of the class so it can be sped up using NUMBA."""

        cpn_times = np.array(cpn_times)
        cpn_flows = np.array(cpn_flows)

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
            face_amount,
            self.sigma,
            self.a,
            self.qq,
            self.pu,
            self.pm,
            self.pd,
            self.r_t,
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

    def df_tree(self, t_mat):
        """Discount factor as seen from now to time t_mat as long as the time
        is on the tree grid."""

        if t_mat == 0.0:
            return 1.0

        _, num_nodes = self.qq.shape
        fn1 = t_mat / self.dt
        fn2 = float(int(t_mat / self.dt))
        if abs(fn1 - fn2) > 1e-6:
            raise FinError("Time not on tree time grid")

        time_step = int(t_mat / self.dt) + 1

        p = 0.0
        for i in range(0, num_nodes):
            ad = self.qq[time_step, i]
            p += ad
        zero_rate = -np.log(p) / t_mat
        return p, zero_rate

    ####################################################################################

    def build_tree(self, tree_mat, df_times, df_values):
        """Build the trinomial tree."""

        if isinstance(df_times, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(df_values, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        # I wish to add on an additional time to the tree so that the second
        # last time corresponds to a maturity tree_mat. For this reason I scale
        # up the maturity date of the tree as follows
        tree_maturity = tree_mat * (self.num_time_steps + 1) / self.num_time_steps

        # The vector of times goes out to this maturity
        tree_times = np.linspace(0.0, tree_maturity, self.num_time_steps + 2)
        self.tree_times = tree_times

        df_tree = np.zeros(shape=self.num_time_steps + 2)
        df_tree[0] = 1.0

        for i in range(1, self.num_time_steps + 2):
            t = tree_times[i]
            df_tree[i] = _uinterpolate(t, df_times, df_values, INTERP_TYPE_VALUE)

        self.df_times = df_times
        self.dfs = df_values

        self.qq, self.pu, self.pm, self.pd, self.r_t, self.dt = build_tree_fast(
            self.a, self.sigma, tree_times, self.num_time_steps, df_tree
        )

        return

    ####################################################################################

    def __repr__(self):
        """Return string with class details."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Sigma", self.sigma)
        s += label_to_string("a", self.a)
        s += label_to_string("num_time_steps", self.num_time_steps)
        s += label_to_string("european_calc_types", self.european_calc_type)
        return s
