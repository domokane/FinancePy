##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum
import numpy as np
from scipy import optimize
from numba import njit
from math import ceil

from ..utils.error import FinError
from ..utils.math import N, accrued_interpolator
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string
from ..utils.global_types import FinExerciseTypes
from ..utils.global_vars import gSmall

interp = InterpTypes.FLAT_FWD_RATES.value

small = 1e-10

###############################################################################
# TODO: Put Jamshidian code into Numba to get speed up
###############################################################################

###############################################################################
# dr = (theta(t) - r) dt + sigma * dW
###############################################################################
###############################################################################


class FinHWEuropeanCalcType(Enum):
    JAMSHIDIAN = 1,
    EXPIRY_ONLY = 2,
    EXPIRY_TREE = 3

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


@njit(fastmath=True, cache=True)
def p_fast(t, T, Rt, delta, pt, ptd, pT, _sigma, _a):
    """ Forward discount factor as seen at some time t which may be in the
    future for payment at time T where Rt is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor. """

    if abs(_a) < small:
        _a = small

    BtT = (1.0 - np.exp(-_a*(T-t)))/_a
    BtDelta = (1.0 - np.exp(-_a * delta))/_a

    term1 = np.log(pT/pt) - (BtT/BtDelta) * np.log(ptd/pt)

    term2 = (_sigma**2)*(1.0-np.exp(-2.0*_a*t)) \
        * BtT * (BtT - BtDelta)/(4.0*_a)

    logAhat = term1 - term2
    BhattT = (BtT/BtDelta) * delta
    p = np.exp(logAhat - BhattT * Rt)
    return p

###############################################################################


@njit(fastmath=True, cache=True)
def build_tree_fast(a, sigma, tree_times, num_time_steps, discount_factors):
    """ Fast tree construction using Numba. """
    treeMaturity = tree_times[-1]
    dt = treeMaturity / (num_time_steps+1)
    dR = sigma * np.sqrt(3.0 * dt)
    jmax = ceil(0.1835/(a * dt))
    N = jmax

    pu = np.zeros(shape=(2*jmax+1))
    pm = np.zeros(shape=(2*jmax+1))
    pd = np.zeros(shape=(2*jmax+1))

    # The short rate goes out one step extra to have the final short rate
    rt = np.zeros(shape=(num_time_steps+2, 2*jmax+1))

    # probabilities start at time 0 and go out to one step before T
    # Branching is simple trinomial out to time step m=1 after which
    # the top node and bottom node connect internally to two lower nodes
    # and two upper nodes respectively. The probabilities only depend on j

    for j in range(-jmax, jmax+1):
        ajdt = a*j*dt
        jN = j + N
        if j == jmax:
            pu[jN] = 7.0/6.0 + 0.50*(ajdt*ajdt - 3.0*ajdt)
            pm[jN] = -1.0/3.0 - ajdt*ajdt + 2.0*ajdt
            pd[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt - ajdt)
        elif j == -jmax:
            pu[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt + ajdt)
            pm[jN] = -1.0/3.0 - ajdt*ajdt - 2.0*ajdt
            pd[jN] = 7.0/6.0 + 0.50*(ajdt*ajdt + 3.0*ajdt)
        else:
            pu[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt - ajdt)
            pm[jN] = 2.0/3.0 - ajdt*ajdt
            pd[jN] = 1.0/6.0 + 0.50*(ajdt*ajdt + ajdt)

    # Arrow-Debreu array
    Q = np.zeros(shape=(num_time_steps+2, 2*N+1))

    # This is the drift adjustment to ensure no arbitrage at each time
    alpha = np.zeros(num_time_steps+1)

    # Time zero is trivial for the Arrow-Debreu price
    Q[0, N] = 1.0

    # Big loop over time steps
    for m in range(0, num_time_steps + 1):

        nm = min(m, jmax)
        sumQZ = 0.0
        for j in range(-nm, nm+1):
            rdt = j*dR*dt
            sumQZ += Q[m, j+N] * np.exp(-rdt)
        alpha[m] = np.log(sumQZ/discount_factors[m+1]) / dt

        for j in range(-nm, nm+1):
            jN = j + N
            rt[m, jN] = alpha[m] + j*dR

        # Loop over all nodes at time m to calculate next values of Q
        for j in range(-nm, nm+1):
            jN = j + N
            rdt = rt[m, jN] * dt
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

###############################################################################


@njit(fastmath=True, cache=True)
def american_bond_option_tree_fast(texp,
                                   strike_price,
                                   face_amount,
                                   coupon_times,
                                   couponAmounts,
                                   exercise_typeInt,
                                   _sigma,
                                   _a,
                                   _Q,
                                   _pu, _pm, _pd,
                                   _rt,
                                   _dt,
                                   _tree_times,
                                   _df_times, _df_values):
    """ Value an option on a bond with coupons that can have European or
    American exercise. Some minor issues to do with handling coupons on
    the option expiry date need to be solved. """

    DEBUG = False
    if DEBUG:
        print("Entering AmerBondOption")
        print("coupon Times", coupon_times)
        print("Coupon Amounts", couponAmounts)

    num_time_steps, num_nodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * dt))
    expiryStep = int(texp/dt + 0.50)

    ###########################################################################

    # Want to add coupons before expiry to the grid so that we can value
    # their impact on the decision to exercise the option early
    treeFlows = np.zeros(num_time_steps)
    num_coupons = len(coupon_times)

    # Flows that fall on the expiry date included. The tree only goes out to
    # the expiry date so coupons after this date do not go onto the tree.
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        if tcpn <= texp:
            n = int(tcpn/dt + 0.50)
            ttree = _tree_times[n]
            df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
            df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
            treeFlows[n] += couponAmounts[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    # I start the tree with the previous coupon time and amount (does not matter)
    mapped_times = np.zeros(0)   # CHANGE
    mapped_amounts = np.zeros(0)  # CHANGE
    for n in range(0, len(_tree_times)):
        if treeFlows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, treeFlows[n])

    # Need future cash flows which are exact time and size for accrued at texp
    for n in range(0, num_coupons):
        if coupon_times[n] > texp:
            mapped_times = np.append(mapped_times, coupon_times[n])
            mapped_amounts = np.append(mapped_amounts, couponAmounts[n])

    if DEBUG:
        print("MAPPED TIMES", mapped_times)
        print("MAPPED AMOUNTS", mapped_amounts)

    ###########################################################################

    accrued = np.zeros(num_time_steps)
    for m in range(0, expiryStep+1):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, coupon_times, couponAmounts)
        accrued[m] *= face_amount

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face_amount

    if DEBUG:
        for i in range(0, expiryStep+1):
            print(i, treeFlows[i], accrued[i])

    ###########################################################################

    call_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    put_option_values = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    ptexp = _uinterpolate(texp, _df_times, _df_values, interp)
    ptdelta = _uinterpolate(texp+dt, _df_times, _df_values, interp)

    cpn = 0.0
    zcb = 0.0

    ###########################################################################
    # As the HW model has a closed form solution for the bond price, I use
    # this fact to calculate the bond price at expiry on the tree nodes
    ###########################################################################

    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + jmax
        rt = _rt[expiryStep, kN]
        bondPrice = 0.0
        for i in range(0, num_coupons):
            tflow = coupon_times[i]
            if tflow >= texp:
                ptflow = _uinterpolate(tflow, _df_times, _df_values, interp)
                zcb = p_fast(texp, tflow, rt, dt, ptexp, ptdelta, ptflow,
                             _sigma, _a)
                cpn = couponAmounts[i]
                bondPrice += cpn * face_amount * zcb

        bondPrice += zcb * face_amount

        # The flow on this date has been added
        bond_values[expiryStep, kN] = bondPrice

    # Now consider exercise of the option on the expiry date
    nm = min(expiryStep, jmax)
    for k in range(-nm, nm+1):
        kN = k + jmax
        full_price = bond_values[expiryStep, kN]
        clean_price = full_price - accrued[expiryStep]
        callExercise = max(clean_price - strike_price, 0.0)
        putExercise = max(strike_price - clean_price, 0.0)
        call_option_values[expiryStep, kN] = callExercise
        put_option_values[expiryStep, kN] = putExercise

    m = expiryStep

    if DEBUG:
        print("-----------------------------------------")
        print("EXP", _tree_times[m], accrued[m], full_price, clean_price,
              callExercise, putExercise)

#        print(kN, bond_values[expiryStep, kN], "CLEAN", clean_price)
#        print("EXPIRY DATE", kN, clean_price, accrued[expiryStep], strike_price)

    # Now step back to today considering exercise at expiry and before
    for m in range(expiryStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face_amount

        for k in range(-nm, nm+1):
            kN = k + jmax
            r = _rt[m, kN]
            df = np.exp(-r*dt)

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
def bermudan_swaption_tree_fast(texp, tmat, strike_price, face_amount,
                                coupon_times, coupon_flows,
                                exercise_typeInt,
                                _df_times, _df_values,
                                _tree_times, _Q, _pu, _pm, _pd, _rt, _dt, _a):
    """ Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to multiple exercise
    times we need to extend tree out to bond maturity and take into account
    cash flows through time. """

    num_time_steps, num_nodes = _Q.shape
    jmax = ceil(0.1835/(_a * _dt))
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    fixed_legFlows = np.zeros(num_time_steps)
    float_leg_values = np.zeros(num_time_steps)
    num_coupons = len(coupon_times)

    # Tree flows go all the way out to the bond maturity date
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        n = int(round(tcpn/_dt, 0))
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        fixed_legFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree
        float_leg_values[n] = strike_price * df_flow / df_tree

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

    ###########################################################################

    # The value of the swap at each time and node. Principal is exchanged.
    fixed_leg_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a payer swap
    pay_values = np.zeros(shape=(num_time_steps, num_nodes))
    # The value of the option to enter into a receiver swap
    rec_values = np.zeros(shape=(num_time_steps, num_nodes))

    # Start with the value of the bond at maturity
    for k in range(0, num_nodes):
        flow = 1.0 + fixed_legFlows[maturityStep]
        fixed_leg_values[maturityStep, k] = flow * face_amount

    N = jmax

    # Now step back to today considering early exercise
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

            elif exercise_typeInt == 2 and flow > gSmall and m >= expiryStep:

                pay_values[m, kN] = max(payExercise, holdPay)
                rec_values[m, kN] = max(recExercise, holdRec)

            elif exercise_typeInt == 3 and m >= expiryStep:

                pay_values[m, kN] = max(payExercise, holdPay)
                rec_values[m, kN] = max(recExercise, holdRec)

                # Need to define floating value on all grid dates

                raise FinError("American optionality not tested.")

    return pay_values[0, jmax], rec_values[0, jmax]

###############################################################################
# TODO: CHECK ACCRUED AND COUPONS TO SEE IF IT WORKS FOR LOW TREE STEPS
###############################################################################


@njit(fastmath=True, cache=True)
def callable_puttable_bond_tree_fast(coupon_times, coupon_flows,
                                     call_times, call_prices,
                                     put_times, put_prices, face,
                                     _sigma, _a, _Q,  # IS SIGMA USED ?
                                     _pu, _pm, _pd, _rt, _dt, _tree_times,
                                     _df_times, _df_values):
    """ Value an option on a bond with coupons that can have European or
    American exercise. Some minor issues to do with handling coupons on
    the option expiry date need to be solved. """

#    print("Coupon Times:", coupon_times)
#    print("Coupon Flows:", coupon_flows)

#    print("DF Times:", _df_times)
#    print("DF Values:", _df_values)

    if np.any(coupon_times < 0.0):
        raise FinError("No coupon times can be before the value date.")

    num_time_steps, num_nodes = _Q.shape
    dt = _dt
    jmax = ceil(0.1835/(_a * dt))
    tmat = coupon_times[-1]
    maturityStep = int(tmat/dt + 0.50)

    ###########################################################################
    # Map coupons onto tree while preserving their present value
    ###########################################################################

    treeFlows = np.zeros(num_time_steps)

    num_coupons = len(coupon_times)
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        n = int(round(tcpn/dt, 0))
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        treeFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree

#    print("Tree flows:", treeFlows)

    ###########################################################################
    # Mapped times stores the mapped times and flows and is used to calculate
    # accrued interest in a consistent manner as using actual flows will
    # result in some convergence noise issues as it is inconsistent
    ###########################################################################

    mapped_times = np.array([0.0])
    mapped_amounts = np.array([0.0])

    for n in range(1, len(_tree_times)):
        if treeFlows[n] > 0.0:
            mapped_times = np.append(mapped_times, _tree_times[n])
            mapped_amounts = np.append(mapped_amounts, treeFlows[n])

    accrued = np.zeros(num_time_steps)
    for m in range(0, num_time_steps):
        ttree = _tree_times[m]
        accrued[m] = accrued_interpolator(ttree, mapped_times, mapped_amounts)
        accrued[m] *= face

        # This is a bit of a hack for when the interpolation does not put the
        # full accrued on flow date. Another scheme may work but so does this
        if treeFlows[m] > 0.0:
            accrued[m] = treeFlows[m] * face

    ###########################################################################
    # map call onto tree - must have no calls at high value
    ###########################################################################

    tree_call_value = np.ones(num_time_steps) * face * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(round(call_time/dt, 0))
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    treePutValue = np.zeros(num_time_steps)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(round(put_time/dt, 0))
        treePutValue[n] = put_prices[i]

    ###########################################################################
    # Value the bond by backward induction starting at bond maturity
    ###########################################################################

    callPutBondValues = np.zeros(shape=(num_time_steps, num_nodes))
    bond_values = np.zeros(shape=(num_time_steps, num_nodes))

    DEBUG = True
    if DEBUG:
        df = 1.0
        px = 0.0
        for i in range(0, maturityStep+1):
            flow = treeFlows[i]
            t = _tree_times[i]
            df = _uinterpolate(t, _df_times, _df_values, interp)

            if flow > gSmall:
                pv = flow * df
                px += pv

        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = min(maturityStep, jmax)
    vcall = tree_call_value[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * face
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(-nm, nm+1):
        kN = k + jmax
        bond_values[m, kN] = (1.0 + treeFlows[m]) * face
        callPutBondValues[m, kN] = value

    # Now step back to today considering early put and call
    for m in range(maturityStep-1, -1, -1):
        nm = min(m, jmax)
        flow = treeFlows[m] * face
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


def fwd_full_bond_price(rt, *args):
    """ Price a coupon bearing bond on the option expiry date and return
    the difference from a strike price. This is used in a root search to
    find the future expiry time short rate that makes the bond price equal
    to the option strike price. It is a key step in the Jamshidian bond
    decomposition approach. The strike is a clean price. """

    self = args[0]
    texp = args[1]
    cpn_times = args[2]
    cpn_amounts = args[3]
    df_times = args[4]
    df_values = args[5]
    strike_price = args[6]
    face = args[7]

    dt = 0.001
    tdelta = texp + dt
    ptexp = _uinterpolate(texp, df_times, df_values, interp)
    ptdelta = _uinterpolate(tdelta, df_times, df_values, interp)

#    print("TEXP", texp, ptexp)

    num_flows = len(cpn_times)
    pv = 0.0

    for i in range(1, num_flows):

        tcpn = cpn_times[i]
        cpn = cpn_amounts[i]

        if tcpn > texp:
            ptcpn = _uinterpolate(tcpn, df_times, df_values, interp)
            zcb = p_fast(texp, tcpn, rt, dt, ptexp, ptdelta, ptcpn,
                         self._sigma, self._a)
            pv = pv + zcb * cpn
#            print("TCPN", tcpn, "ZCB", zcb, "CPN", cpn, "PV", pv)

    if tcpn >= texp:
        pv = pv + zcb

#    print("TCPN", tcpn, "ZCB", zcb, "PRI", 1.0, "PV", pv)

    accd = accrued_interpolator(texp, cpn_times, cpn_amounts)
#    print("Accrued:", accd)
#    print("texp:", texp)
#    print("cpn_times:", cpn_times)
#    print("cpn_amounts:", cpn_amounts)

    pv_clean = pv - accd
    obj = face * pv_clean - strike_price

#    print("FWD PRICE", rt, pv, accd, strike_price, obj)
    return obj

###############################################################################


class HWTree():

    def __init__(self,
                 sigma,
                 a,
                 num_time_steps=100,
                 europeanCalcType=FinHWEuropeanCalcType.EXPIRY_TREE):
        """ Constructs the Hull-White rate model. The speed of mean reversion
        a and volatility are passed in. The short rate process is given by
        dr = (theta(t) - ar) * dt  + sigma * dW. The model will switch to use
        Jamshidian's approach where possible unless the useJamshidian flag is
        set to false in which case it uses the trinomial Tree. """

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        if a < 0.0:
            raise FinError("Mean reversion speed parameter should be >= 0.")

        self._sigma = sigma
        self._a = a
        self._num_time_steps = num_time_steps
        self._europeanCalcType = europeanCalcType

        self._Q = None
        self._r = None
        self._tree_times = None
        self._pu = None
        self._pm = None
        self._pd = None
        self._discount_curve = None
        self._treeBuilt = False

###############################################################################

    def option_on_zcb(self,
                      texp, tmat,
                      strike, face_amount,
                      df_times, df_values):
        """ Price an option on a zero coupon bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date. """

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        ptexp = _uinterpolate(texp, df_times, df_values, interp)
        ptmat = _uinterpolate(tmat, df_times, df_values, interp)

        sigma = self._sigma
        a = self._a

        if abs(a) < small:
            a = small

        sigmap = (sigma/a) * (1.0 - np.exp(-a*(tmat-texp)))
        sigmap *= np.sqrt((1.0-np.exp(-2.0*a*texp))/2.0/a)

        if abs(sigmap) < small:
            sigmap = small

        h = np.log((face_amount*ptmat)/(strike * ptexp)) / \
            sigmap + sigmap / 2.0
        callValue = face_amount * ptmat * N(h) - strike * ptexp * N(h - sigmap)
        putValue = strike * ptexp * \
            N(-h + sigmap) - face_amount * ptmat * N(-h)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def european_bond_option_jamshidian(self,
                                        texp,
                                        strike_price,
                                        face,
                                        cpn_times,
                                        cpn_amounts,
                                        df_times,
                                        df_values):
        """ Valuation of a European bond option using the Jamshidian
        deconstruction of the bond into a strip of zero coupon bonds with the
        short rate that would make the bond option be at the money forward. """

#        print(df_times)
#        print(df_values)

        num_coupons = len(cpn_times)

        argtuple = (self, texp, cpn_times, cpn_amounts,
                    df_times, df_values, strike_price, face)

        # Can I improve on this initial guess ?
        x0 = 0.05

        rstar = optimize.newton(fwd_full_bond_price, x0=x0, fprime=None,
                                args=argtuple, tol=1e-10, maxiter=50,
                                fprime2=None)

        # Now we price a series of zero coupon bonds using this short rate
        dt = 1e-6

        ptexp = _uinterpolate(texp, df_times, df_values, interp)
        ptdelta = _uinterpolate(texp+dt, df_times, df_values, interp)

        callValue = 0.0
        putValue = 0.0

        # Adjust strike to handle
        for i in range(0, num_coupons):

            tcpn = cpn_times[i]
            cpn = cpn_amounts[i]

            if tcpn >= texp:  # coupons on the expiry date are included

                ptcpn = _uinterpolate(tcpn, df_times, df_values, interp)

                strike = p_fast(texp, tcpn, rstar, dt, ptexp, ptdelta,
                                ptcpn, self._sigma, self._a)

                v = self.option_on_zcb(texp, tcpn, strike, 1.0,
                                       df_times, df_values)

                call = v['call']
                put = v['put']

                callValue += call * cpn * face
                putValue += put * cpn * face

        callValue += call * face
        putValue += put * face

        return {'call': callValue, 'put': putValue}

###############################################################################

    def european_bond_option_expiry_only(self,
                                         texp,
                                         strike_price,
                                         face_amount,
                                         cpn_times,
                                         cpn_amounts):
        """ Price a European option on a coupon-paying bond using a tree to
        generate short rates at the expiry date and then to use the analytical
        solution of zero coupon bond prices in the HW model to calculate the
        corresponding bond price. User provides bond object and option details.
        """

        dt = self._dt
        tdelta = texp + dt

        ptexp = _uinterpolate(texp, self._df_times, self._dfs, interp)
        ptdelta = _uinterpolate(tdelta, self._df_times, self._dfs, interp)

        _, num_nodes = self._Q.shape
        expiryStep = int(texp/dt+0.50)

        callValue = 0.0
        putValue = 0.0
        num_coupons = len(cpn_times)

        #######################################################################

        for k in range(0, num_nodes):

            q = self._Q[expiryStep, k]
            rt = self._rt[expiryStep, k]

            pv = 0.0

            for i in range(0, num_coupons):

                tcpn = cpn_times[i]
                cpn = cpn_amounts[i]

                if tcpn >= texp:

                    ptcpn = _uinterpolate(tcpn, self._df_times, self._dfs,
                                          interp)

                    zcb = p_fast(texp, tcpn, rt, dt, ptexp, ptdelta, ptcpn,
                                 self._sigma, self._a)

                    pv += cpn * zcb

            pv += zcb

#            print(texp)
#            print(cpn_times)
#            print(cpn_amounts)

            accrued = accrued_interpolator(texp, cpn_times, cpn_amounts)

            pv = pv - accrued
#            print(accrued)

            putPayoff = max(strike_price - pv * face_amount, 0.0)
            callPayoff = max(pv * face_amount - strike_price, 0.0)

            putValue += q * putPayoff
            callValue += q * callPayoff

        #######################################################################

        return {'call': callValue, 'put': putValue}

###############################################################################

    def option_on_zero_coupon_bond_tree(self,
                                        texp,
                                        tmat,
                                        strike_price,
                                        face_amount):
        """ Price an option on a zero coupon bond using a HW trinomial
        tree. The discount curve was already supplied to the tree build. """

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        if self._tree_times is None:
            raise FinError("Tree has not been constructed.")

        if self._tree_times[-1] < texp:
            raise FinError("Tree expiry must be >= option expiry date.")

        dt = self._dt
        tdelta = texp + dt

        ptexp = _uinterpolate(texp, self._df_times, self._dfs, interp)
        ptdelta = _uinterpolate(tdelta, self._df_times, self._dfs, interp)
        ptmat = _uinterpolate(tmat, self._df_times, self._dfs, interp)

        _, num_nodes = self._Q.shape
        expiryStep = int(texp/dt+0.50)

        callValue = 0.0
        putValue = 0.0

        for k in range(0, num_nodes):

            q = self._Q[expiryStep, k]
            rt = self._rt[expiryStep, k]

            zcb = p_fast(texp, tmat,
                         rt, dt, ptexp, ptdelta, ptmat,
                         self._sigma, self._a)

            putPayoff = max(strike_price - zcb * face_amount, 0.0)
            callPayoff = max(zcb * face_amount - strike_price, 0.0)
            putValue += q * putPayoff
            callValue += q * callPayoff

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudan_swaption(self, texp, tmat, strike, face,
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
            = bermudan_swaption_tree_fast(texp, tmat, strike, face,
                                          coupon_times, coupon_flows,
                                          exercise_typeInt,
                                          self._df_times, self._dfs,
                                          self._tree_times, self._Q,
                                          self._pu, self._pm, self._pd,
                                          self._rt,
                                          self._dt, self._a)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def bond_option(self, texp, strike_price, face_amount,
                    coupon_times, coupon_flows, exercise_type):
        """ Value a bond option that can have European or American exercise.
        This is done using a trinomial tree that we extend out to bond 
        maturity. For European bond options, Jamshidian's model is
        faster and is used instead i.e. not this function. """

        exercise_typeInt = option_exercise_types_to_int(exercise_type)

        if exercise_typeInt == 1:

            if self._europeanCalcType == FinHWEuropeanCalcType.JAMSHIDIAN:

                v = self.european_bond_option_jamshidian(texp,
                                                         strike_price,
                                                         face_amount,
                                                         coupon_times,
                                                         coupon_flows,
                                                         self._df_times,
                                                         self._dfs)

                callValue = v['call']
                putValue = v['put']

            elif self._europeanCalcType == FinHWEuropeanCalcType.EXPIRY_ONLY:

                v = self.european_bond_option_expiry_only(texp,
                                                          strike_price,
                                                          face_amount,
                                                          coupon_times,
                                                          coupon_flows)

                callValue = v['call']
                putValue = v['put']

            elif self._europeanCalcType == FinHWEuropeanCalcType.EXPIRY_TREE:

                callValue, putValue \
                    = american_bond_option_tree_fast(texp,
                                                     strike_price, face_amount,
                                                     coupon_times, coupon_flows,
                                                     exercise_typeInt,
                                                     self._sigma, self._a,
                                                     self._Q,
                                                     self._pu, self._pm, self._pd,
                                                     self._rt, self._dt,
                                                     self._tree_times,
                                                     self._df_times, self._dfs)

            else:
                raise FinError("Unknown HW model implementation choice.")

        else:

            callValue, putValue \
                = american_bond_option_tree_fast(texp,
                                                 strike_price, face_amount,
                                                 coupon_times, coupon_flows,
                                                 exercise_typeInt,
                                                 self._sigma, self._a,
                                                 self._Q,
                                                 self._pu, self._pm, self._pd,
                                                 self._rt, self._dt,
                                                 self._tree_times,
                                                 self._df_times, self._dfs)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def callable_puttable_bond_tree(self,
                                    coupon_times,
                                    coupon_flows,
                                    call_times,
                                    call_prices,
                                    put_times,
                                    put_prices,
                                    face_amount):
        """ Value an option on a bond with coupons that can have European or
        American exercise. Some minor issues to do with handling coupons on
        the option expiry date need to be solved. Also this function should be
        moved out of the class so it can be sped up using NUMBA. """

        coupon_times = np.array(coupon_times)
        coupon_flows = np.array(coupon_flows)

        call_times = np.array(call_times)
        put_times = np.array(put_times)

        call_prices = np.array(call_prices)
        put_prices = np.array(put_prices)

        v = callable_puttable_bond_tree_fast(coupon_times, coupon_flows,
                                             call_times, call_prices,
                                             put_times, put_prices,
                                             face_amount,
                                             self._sigma, self._a,
                                             self._Q,
                                             self._pu, self._pm, self._pd,
                                             self._rt, self._dt,
                                             self._tree_times,
                                             self._df_times, self._dfs)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def df_tree(self, tmat):
        """ Discount factor as seen from now to time tmat as long as the time
        is on the tree grid. """

        if tmat == 0.0:
            return 1.0

        _, num_nodes = self._Q.shape
        fn1 = tmat/self._dt
        fn2 = float(int(tmat/self._dt))
        if abs(fn1 - fn2) > 1e-6:
            raise FinError("Time not on tree time grid")

        timeStep = int(tmat / self._dt) + 1

        p = 0.0
        for i in range(0, num_nodes):
            ad = self._Q[timeStep, i]
            p += ad
        zero_rate = -np.log(p)/tmat
        return p, zero_rate

###############################################################################

    def build_tree(self, treeMat, df_times, df_values):
        """ Build the trinomial tree. """

        if isinstance(df_times, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(df_values, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        # I wish to add on an additional time to the tree so that the second
        # last time corresponds to a maturity treeMat. For this reason I scale
        # up the maturity date of the tree as follows
        treeMaturity = treeMat * (self._num_time_steps+1)/self._num_time_steps

        # The vector of times goes out to this maturity
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

    def __repr__(self):
        """ Return string with class details. """

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Sigma", self._sigma)
        s += label_to_string("a", self._a)
        s += label_to_string("num_time_steps", self._num_time_steps)
        s += label_to_string("EuropeanCalcTypes", self._europeanCalcType)
        return s

###############################################################################
