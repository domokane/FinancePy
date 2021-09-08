##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64, int64

from ..utils.error import FinError
from ..utils.math import accrued_interpolator
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.helpers import label_to_string
from ..utils.global_types import FinExerciseTypes
from ..utils.global_vars import gSmall

interp = InterpTypes.FLAT_FWD_RATES.value

###############################################################################
# TODO: PUT CALL PARITY IS NOT EXACTLY OBSERVED FOR BERMUDAN SWAPTIONS WHEN
#       VOL IS TURNED UP. SMALL EFFECT. $3 OUT OF $1m. SEE TESTFINBERMUDANSWAPTIONS
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


@njit(float64(float64, int64, float64[:, :], float64[:, :],
              float64, float64, float64),
      fastmath=True, cache=True)
def f(x0, m, Q, rt, dfEnd, dt, sigma):

    # x is the middle value on the short-rate on the tree
    midm = int(m/2)
    rt[m, midm] = x0

    for i in range(midm, 0, -1):
        rt[m, i-1] = rt[m, i] * np.exp(-2.0 * sigma * np.sqrt(dt))

    for i in range(midm + 1, m + 1, 1):
        rt[m, i] = rt[m, i-1] * np.exp(2.0 * sigma * np.sqrt(dt))

    sumInner = 0.0
    for i in range(0, m+1):
        r = rt[m, i]
        nextPeriodDf = (1.0 / ((1.0 + r)**dt))
        q = Q[m, i]
        sumInner += q * nextPeriodDf

    obj_fn = sumInner - dfEnd
    return obj_fn

###############################################################################


@njit(float64(float64, int64, float64[:, :], float64[:, :], float64, float64,
              float64), fastmath=True, cache=True)
def search_root(x0, m, Q, rt, dfEnd, dt, sigma):

    #    print("Searching for root", x0)
    max_iter = 10
    max_error = 1e-8

    x1 = x0 * 1.0001
    f0 = f(x0, m, Q, rt, dfEnd, dt, sigma)
    f1 = f(x1, m, Q, rt, dfEnd, dt, sigma)

    for _ in range(0, max_iter):

        df = f1 - f0

        if df == 0.0:
            raise FinError("Search for alpha fails due to zero derivative")

        x = x1 - f1 * (x1-x0)/df
        x0, f0 = x1, f1
        x1 = x
        f1 = f(x1, m, Q, rt, dfEnd, dt, sigma)

        if (abs(f1) <= max_error):
            return x1

    raise FinError("Search root deriv FAILED to find alpha.")

###############################################################################


@njit(fastmath=True, cache=True)
def bermudan_swaption_tree_fast(texp, tmat,
                                strike_price, face_amount,
                                coupon_times, coupon_flows,
                                exercise_typeInt,
                                _df_times, _df_values,
                                _tree_times,
                                _Q, _rt, _dt):
    """ Option to enter into a swap that can be exercised on coupon payment
    dates after the start of the exercise period. Due to non-analytical bond
    price we need to extend tree out to bond maturity and take into account
    cash flows through time. """

    pu = 0.50
    pd = 0.50

    ###########################################################################

    num_time_steps, num_nodes = _Q.shape
    expiryStep = int(texp/_dt + 0.50)
    maturityStep = int(tmat/_dt + 0.50)

    ###########################################################################

    fixed_legFlows = np.zeros(num_time_steps)
    float_leg_values = np.zeros(num_time_steps)
    num_coupons = len(coupon_times)

    # Tree flows go all the way out to the swap maturity date
    for i in range(0, num_coupons):
        tcpn = coupon_times[i]
        n = int(tcpn/_dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        fixed_legFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree
        float_leg_values[n] = strike_price  # * df_flow / df_tree

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

    # Now step back to today considering early exercise on coupon dates
    for m in range(maturityStep-1, -1, -1):
        nm = m
        flow = fixed_legFlows[m] * face_amount

        for k in range(0, nm+1):
            rt = _rt[m, k]
            df = np.exp(- rt * _dt)

            vu = fixed_leg_values[m+1, k+1]
            vd = fixed_leg_values[m+1, k]
            v = (pu*vu + pd*vd) * df

            fixed_leg_values[m, k] = v
            fixed_leg_values[m, k] += flow

            vpay = 0.0
            vrec = 0.0

            vu = pay_values[m+1, k+1]
            vd = pay_values[m+1, k]
            vpay = (pu*vu + pd*vd) * df
            pay_values[m, k] = vpay

            vu = rec_values[m+1, k+1]
            vd = rec_values[m+1, k]
            vrec = (pu*vu + pd*vd) * df
            rec_values[m, k] = vrec

            holdPay = pay_values[m, k]
            holdRec = rec_values[m, k]

            # The floating value is clean and so must be the fixed value
            fixed_leg_value = fixed_leg_values[m, k] - accrued[m]
            float_leg_value = float_leg_values[m]

            payExercise = max(float_leg_value - fixed_leg_value, 0.0)
            recExercise = max(fixed_leg_value - float_leg_value, 0.0)

            if m == expiryStep:

                pay_values[m, k] = max(payExercise, holdPay)
                rec_values[m, k] = max(recExercise, holdRec)

            elif exercise_typeInt == 2 and flow > gSmall and m > expiryStep:

                pay_values[m, k] = max(payExercise, holdPay)
                rec_values[m, k] = max(recExercise, holdRec)

            elif exercise_typeInt == 3 and m > expiryStep:

                raise FinError("American optionality not completed.")

                # Need to define floating value on all grid dates

                pay_values[m, k] = max(payExercise, holdPay)
                rec_values[m, k] = max(recExercise, holdRec)

    return pay_values[0, 0], rec_values[0, 0]

###############################################################################


@njit(fastmath=True, cache=True)
def american_bond_option_tree_fast(texp, tmat,
                                   strike_price, face_amount,
                                   coupon_times, coupon_flows,
                                   exercise_typeInt,
                                   _df_times, _df_values,
                                   _tree_times, _Q,
                                   _rt, _dt):
    """ Option to buy or sell bond at a specified strike price that can be
    exercised over the exercise period depending on choice of exercise type.
    Due to non-analytical bond price we need to extend tree out to bond
    maturity and take into account cash flows through time. """

    DEBUG = False

    ###########################################################################

    pu = 0.50
    pd = 0.50

    if DEBUG:
        print("COUPON TIMES", coupon_times)
        print("COUPON AMOUNTS", coupon_flows)

    ###########################################################################

    num_time_steps, num_nodes = _Q.shape
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
            print(coupon_times)
            raise FinError("Coupon times must be positive.")

        n = int(tcpn/_dt + 0.50)
        ttree = _tree_times[n]
        df_flow = _uinterpolate(tcpn, _df_times, _df_values, interp)
        df_tree = _uinterpolate(ttree, _df_times, _df_values, interp)
        treeFlows[n] += coupon_flows[i] * 1.0 * df_flow / df_tree

    ###########################################################################
    # mapped_times = np.zeros(0)   # CHANGE
    # mapped_amounts = np.zeros(0)  # CHANGE
    # for n in range(0, len(_tree_times)):
    #     if treeFlows[n] > 0.0:
    #         mapped_times = np.append(mapped_times, _tree_times[n])
    #         mapped_amounts = np.append(mapped_amounts, treeFlows[n])
    # if DEBUG:
    #     print("MAPPED TIMES", mapped_times)
    #     print("MAPPED AMOUNTS", mapped_amounts)
    #
    #    if mapped_times[0] > gSmall:
    #        raise FinError("Mapped times [0] must be <= 0 for first coupon > 0")
    #
    ###########################################################################

    if DEBUG:
        print(_tree_times)
        print(treeFlows)

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

    ###########################################################################

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

        nm = m
        flow = treeFlows[m] * face_amount

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r * _dt)

            vu = bond_values[m+1, k+1]
            vd = bond_values[m+1, k]
            v = (pu*vu + pd*vd) * df

            bond_values[m, k] = v
            bond_values[m, k] += flow
            full_price = bond_values[m, k]
            clean_price = full_price - accrued[m]

        if DEBUG:
            print(m, _tree_times[m], accrued[m], full_price, clean_price, 0, 0)

    # Now step back to today from the expiry date considering early exercise
    for m in range(expiryStep, -1, -1):
        nm = m
        flow = treeFlows[m] * face_amount

        for k in range(0, nm+1):

            r = _rt[m, k]
            df = np.exp(-r*_dt)

            vu = bond_values[m+1, k+1]
            vd = bond_values[m+1, k]
            v = (pu*vu + pd*vd) * df

            bond_values[m, k] = v
            bond_values[m, k] += flow

            vcall = 0.0
            vput = 0.0

            vu = call_option_values[m+1, k+1]
            vd = call_option_values[m+1, k]
            vcall = (pu*vu + pd*vd) * df

            call_option_values[m, k] = vcall

            vu = put_option_values[m+1, k+1]
            vd = put_option_values[m+1, k]
            vput = (pu*vu + pd*vd) * df

            put_option_values[m, k] = vput

            full_price = bond_values[m, k]
            clean_price = full_price - accrued[m]
            callExercise = max(clean_price - strike_price, 0.0)
            putExercise = max(strike_price - clean_price, 0.0)

            holdCall = call_option_values[m, k]
            holdPut = put_option_values[m, k]

            if m == expiryStep:

                call_option_values[m, k] = max(callExercise, holdCall)
                put_option_values[m, k] = max(putExercise, holdPut)

            elif exercise_typeInt == 3 and m < expiryStep:

                call_option_values[m, k] = max(callExercise, holdCall)
                put_option_values[m, k] = max(putExercise, holdPut)

        if DEBUG:
            print(m, _tree_times[m], accrued[m], full_price, clean_price,
                  callExercise, putExercise)

    return call_option_values[0, 0], put_option_values[0, 0]

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

    pu = 0.50
    pd = 0.50

    #######################################################################
    num_time_steps, num_nodes = _Q.shape
    dt = _dt
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
            px += flow * df
        px += df

    ###########################################################################
    # Now step back to today considering early exercise
    ###########################################################################

    m = maturityStep
    nm = maturityStep
    vcall = tree_call_value[m]
    vput = treePutValue[m]
    vhold = (1.0 + treeFlows[m]) * face_amount
    vclean = vhold - accrued[m]
    value = min(max(vclean, vput), vcall) + accrued[m]

    for k in range(0, nm+1):
        bond_values[m, k] = (1.0 + treeFlows[m]) * face_amount
        callPutBondValues[m, k] = value

    for m in range(maturityStep-1, -1, -1):
        nm = m
        flow = treeFlows[m] * face_amount
        vcall = tree_call_value[m]
        vput = treePutValue[m]

        for k in range(0, nm+1):

            rt = _rt[m, k]
            df = np.exp(-rt*dt)

            vu = bond_values[m+1, k+1]
            vd = bond_values[m+1, k]
            v = (pu*vu + pd*vd) * df

            bond_values[m, k] = v
            bond_values[m, k] += flow

            vu = callPutBondValues[m+1, k+1]
            vd = callPutBondValues[m+1, k]

            vhold = (pu*vu + pd*vd) * df
            # Need to make add on coupons paid if we hold
            vhold = vhold + flow
            value = min(max(vhold - accrued[m], vput), vcall) + accrued[m]
            callPutBondValues[m, k] = value

    return {'bondwithoption': callPutBondValues[0, 0],
            'bondpure': bond_values[0, 0]}

###############################################################################
###############################################################################


@njit(cache=True, fastmath=True)
def build_tree_fast(sigma, tree_times, num_time_steps, discount_factors):
    # Unlike the BK and HW Trinomial trees, this Tree is packed into the lower
    # diagonal of a square matrix because of its binomial nature. This means
    # that the indexing of the arrays is different.

    treeMaturity = tree_times[-1]
    dt = treeMaturity / (num_time_steps+1)

    # The short rate goes out one step extra to have the final short rate
    # as it follows HW code but I am not sure this is needed. EXAMINE
    # NOTE HW code uses this to have short rate at expiry so it can use
    # analytical solutions for the zero coupon bond price
    # This is the BDT model so x = log(r)

    # Some implementations use continuous compounding. Does not affect results.
    # Useful if you want to make direct comparisons of the tree.
    CONT_COMPOUNDED = True

    Q = np.zeros(shape=(num_time_steps+2, num_time_steps+2))
    rt = np.zeros(shape=(num_time_steps+2, num_time_steps+2))

    if CONT_COMPOUNDED:
        r0 = -np.log(discount_factors[1])/dt
    else:
        r0 = (1.0 / discount_factors[1] - 1.0)/dt

    rt[0, 0] = r0
    Q[0, 0] = 1.0

    if CONT_COMPOUNDED:
        Q[1, 0] = 0.50 * np.exp(-rt[0, 0] * dt)
        Q[1, 1] = 0.50 * np.exp(-rt[0, 0] * dt)
    else:
        Q[1, 0] = 0.50 / ((1.0 + rt[0, 0]) ** dt)
        Q[1, 1] = 0.50 / ((1.0 + rt[0, 0]) ** dt)

    # The short rate goes out one step extra to have the final short rate
    # as it follows HW code but I am not sure this is needed. EXAMINE
    for m in range(1, num_time_steps+1):

        dfEnd = discount_factors[m+1]
        search_root(r0, m, Q, rt, dfEnd, dt, sigma)

        if CONT_COMPOUNDED:
            Q[m+1, 0] = 0.50 * Q[m, 0] * np.exp(-rt[m, 0] * dt)
        else:
            Q[m+1, 0] = 0.50 * Q[m, 0] / ((1.0 + rt[m, 0]) ** dt)

        for n in range(1, m+1):

            if CONT_COMPOUNDED:
                Q[m+1, n] = 0.50 * Q[m, n-1] * np.exp(-rt[m, n-1] * dt) \
                    + 0.50 * Q[m, n] * np.exp(-rt[m, n] * dt)
            else:
                Q[m+1, n] = 0.50 * Q[m, n-1] / ((1.0 + rt[m, n-1])**dt) \
                    + 0.50 * Q[m, n] / ((1.0 + rt[m, n])**dt)

        if CONT_COMPOUNDED:
            Q[m+1, m+1] = 0.50 * Q[m, m] * np.exp(-rt[m, m] * dt)
        else:
            Q[m+1, m+1] = 0.50 * Q[m, m] / ((1.0 + rt[m, m])**dt)

    return (Q, rt, dt)

###############################################################################


class BDTTree():

    def __init__(self,
                 sigma: float,
                 num_time_steps: int = 100):
        """ Constructs the Black-Derman-Toy rate model in the case when the
        volatility is assumed to be constant. The short rate process simplifies
        and is given by d(log(r)) = theta(t) * dt + sigma * dW. Althopugh """

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        self._sigma = sigma

        if num_time_steps < 3:
            raise FinError("Drift fitting requires at least 3 time steps.")

        self._num_time_steps = num_time_steps

        self._Q = None
        self._rt = None
        self._tree_times = None
        self._pu = 0.50
        self._pd = 0.50
        self._discount_curve = None

###############################################################################

    def build_tree(self, treeMat, df_times, df_values):

        if isinstance(df_times, np.ndarray) is False:
            raise FinError("DF TIMES must be a numpy vector")

        if isinstance(df_values, np.ndarray) is False:
            raise FinError("DF VALUES must be a numpy vector")

        interp = InterpTypes.FLAT_FWD_RATES.value

        treeMaturity = treeMat * (self._num_time_steps+1)/self._num_time_steps
        tree_times = np.linspace(0.0, treeMaturity, self._num_time_steps + 2)
        self._tree_times = tree_times

        dfTree = np.zeros(shape=(self._num_time_steps+2))
        dfTree[0] = 1.0

        for i in range(1, self._num_time_steps+2):
            t = tree_times[i]
            dfTree[i] = _uinterpolate(t, df_times, df_values, interp)

        self._df_times = df_times
        self._dfs = df_values

        self._Q, self._rt, self._dt \
            = build_tree_fast(self._sigma,
                              tree_times, self._num_time_steps, dfTree)

        return

###############################################################################

    def bond_option(self, texp, strike_price, face_amount,
                    coupon_times, coupon_flows, exercise_type):
        """ Value a bond option that can have European or American exercise
        using the Black-Derman-Toy model. The model uses a binomial tree. """

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
                                             self._rt,
                                             self._dt)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def bermudan_swaption(self, texp, tmat, strike, face_amount,
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
                                          strike, face_amount,
                                          coupon_times, coupon_flows,
                                          exercise_typeInt,
                                          self._df_times, self._dfs,
                                          self._tree_times, self._Q,
                                          self._rt,
                                          self._dt)

        return {'pay': payValue, 'rec': recValue}

###############################################################################

    def callable_puttable_bond_tree(self,
                                    coupon_times, coupon_flows,
                                    call_times, call_prices,
                                    put_times, put_prices,
                                    face_amount):
        """ Option that can be exercised at any time over the exercise period.
        Due to non-analytical bond price we need to extend tree out to bond
        maturity and take into account cash flows through time. """

        call_times = np.array(call_times)
        put_times = np.array(put_times)

        call_prices = np.array(call_prices)
        put_prices = np.array(put_prices)

        v = callable_puttable_bond_tree_fast(coupon_times, coupon_flows,
                                             call_times, call_prices,
                                             put_times, put_prices, face_amount,
                                             self._sigma,
                                             self._Q,
                                             self._rt, self._dt,
                                             self._tree_times,
                                             self._df_times, self._dfs)

        return {'bondwithoption': v['bondwithoption'],
                'bondpure': v['bondpure']}

###############################################################################

    def __repr__(self):
        """ Return string with class details. """

        s = "Black-Derman-Toy Model\n"
        s += label_to_string("Sigma", self._sigma)
        s += label_to_string("num_time_steps", self._num_time_steps)
        return s

###############################################################################
