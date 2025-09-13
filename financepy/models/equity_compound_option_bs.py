from math import exp, log, sqrt

import numpy as np
from scipy import optimize
from numba import njit

from ..utils.math import normcdf, phi2
from ..utils.global_vars import G_SMALL
from ..utils.error import FinError
from ..utils.global_types import OptionTypes

from .black_scholes_analytic import bs_value

########################################################################################


def _f(s0, *args) -> float:
    opt2_type = args[0]
    t2: float = args[1]
    k2: float = args[2]
    t1: float = args[3]
    r: float = args[4]
    q: float = args[5]
    vol: float = args[6]
    k1: float = args[7]  # Target

    if s0 <= 0.0:
        raise FinError("Unable to solve for stock price that fits k_1")

    tau: float = t2 - t1
    opt_value: float = bs_value(s0, tau, k2, r, q, vol, opt2_type.value)
    obj_fn: float = opt_value - k1

    return obj_fn


########################################################################################



from typing import Any

@njit(fastmath=True, cache=True)
def value_cmpd_once(
    s: float,
    r: float,
    q: float,
    volatility: float,
    t1: float,
    t2: float,
    opt_type1: Any,  # OptionTypes, but Numba doesn't support Enums in signatures
    opt_type2: Any,  # OptionTypes
    k1: float,
    k2: float,
    num_steps: int
) -> np.ndarray:

    if num_steps < 3:
        num_steps = 3

    # TODO EUROPEAN call-put works but AMERICAN call-put needs to be tested

    # Need equally spaced time intervals for a recombining tree
    # Downside is that we may not measure periods exactly
    dt = t2 / num_steps
    num_steps1 = int(t1 / dt)
    num_steps2 = num_steps - num_steps1
    dt1 = dt
    dt2 = dt

    # print("T1:",t1,"T2:",t2,"dt:",dt,"N1*dt",num_steps1*dt,"N*dt",num_steps*dt)

    # the number of nodes on the tree
    num_nodes = (num_steps + 1) * (num_steps + 2) // 2
    option_vals = np.zeros(num_nodes)

    u1 = np.exp(volatility * np.sqrt(dt))
    d1 = 1.0 / u1
    u2 = np.exp(volatility * np.sqrt(dt))
    d2 = 1.0 / u2

    probs = np.zeros(num_steps)
    period_dfs = np.zeros(num_steps)

    # store time independent information for later use in tree
    for i_time in range(0, num_steps1):
        a1 = np.exp((r - q) * dt1)
        probs[i_time] = (a1 - d1) / (u1 - d1)
        period_dfs[i_time] = np.exp(-r * dt1)

    for i_time in range(num_steps1, num_steps):
        a2 = np.exp((r - q) * dt2)
        probs[i_time] = (a2 - d2) / (u2 - d2)
        period_dfs[i_time] = np.exp(-r * dt2)

    stock_vals = np.zeros(num_nodes)
    stock_vals[0] = s
    s_low = s

    for i_time in range(1, num_steps1 + 1):
        s_low *= d1
        s = s_low
        for i_node in range(0, i_time + 1):
            index = i_time * (i_time + 1) // 2
            stock_vals[index + i_node] = s
            s = s * (u1 * u1)

    for i_time in range(num_steps1 + 1, num_steps + 1):
        s_low *= d2
        s = s_low
        for i_node in range(0, i_time + 1):
            index = i_time * (i_time + 1) // 2
            stock_vals[index + i_node] = s
            s = s * (u2 * u2)

    # work backwards by first setting values at expiry date t2
    index = num_steps * (num_steps + 1) // 2

    for i_node in range(0, i_time + 1):
        s = stock_vals[index + i_node]
        if (
            opt_type2 == OptionTypes.EUROPEAN_CALL
            or opt_type2 == OptionTypes.AMERICAN_CALL
        ):
            option_vals[index + i_node] = max(s - k2, 0.0)
        elif (
            opt_type2 == OptionTypes.EUROPEAN_PUT
            or opt_type2 == OptionTypes.AMERICAN_PUT
        ):
            option_vals[index + i_node] = max(k2 - s, 0.0)

    # begin backward steps from expiry at t2 to first expiry at time t1
    for i_time in range(num_steps - 1, num_steps1, -1):
        index = i_time * (i_time + 1) // 2
        for i_node in range(0, i_time + 1):
            s = stock_vals[index + i_node]
            next_index = (i_time + 1) * (i_time + 2) // 2
            next_node_dn = next_index + i_node
            next_node_up = next_index + i_node + 1
            v_up = option_vals[next_node_up]
            v_dn = option_vals[next_node_dn]
            future_exp_val = probs[i_time] * v_up
            future_exp_val += (1.0 - probs[i_time]) * v_dn
            hold_value = period_dfs[i_time] * future_exp_val

            exercise_value = 0.0  # NUMBA NEEDS HELP TO DETERMINE THE TYPE

            if opt_type1 == OptionTypes.AMERICAN_CALL:
                exercise_value = max(s - k2, 0.0)
            elif opt_type1 == OptionTypes.AMERICAN_PUT:
                exercise_value = max(k2 - s, 0.0)

            option_vals[index + i_node] = max(exercise_value, hold_value)

    # Now do payoff at the end of the first expiry period at t1
    i_time = num_steps1
    index = i_time * (i_time + 1) // 2

    for i_node in range(0, i_time + 1):
        s = stock_vals[index + i_node]
        next_index = (i_time + 1) * (i_time + 2) // 2
        next_node_dn = next_index + i_node
        next_node_up = next_index + i_node + 1
        v_up = option_vals[next_node_up]
        v_dn = option_vals[next_node_dn]
        future_exp_val = probs[i_time] * v_up
        future_exp_val += (1.0 - probs[i_time]) * v_dn
        hold_value = period_dfs[i_time] * future_exp_val

        if (
            opt_type1 == OptionTypes.EUROPEAN_CALL
            or opt_type1 == OptionTypes.AMERICAN_CALL
        ):
            option_vals[index + i_node] = max(hold_value - k1, 0.0)
        elif (
            opt_type1 == OptionTypes.EUROPEAN_PUT
            or opt_type1 == OptionTypes.AMERICAN_PUT
        ):
            option_vals[index + i_node] = max(k1 - hold_value, 0.0)

    # begin backward steps from t1 expiry to value date
    for i_time in range(num_steps1 - 1, -1, -1):

        index = i_time * (i_time + 1) // 2

        for i_node in range(0, i_time + 1):

            s = stock_vals[index + i_node]

            next_index = (i_time + 1) * (i_time + 2) // 2

            next_node_dn = next_index + i_node
            next_node_up = next_index + i_node + 1

            v_up = option_vals[next_node_up]
            v_dn = option_vals[next_node_dn]

            future_exp_val = probs[i_time] * v_up
            future_exp_val += (1.0 - probs[i_time]) * v_dn

            hold_value = period_dfs[i_time] * future_exp_val

            exercise_value = 0.0  # NUMBA NEEDS HELP TO DETERMINE THE TYPE

            if opt_type1 == OptionTypes.AMERICAN_CALL:
                exercise_value = max(hold_value - k1, 0.0)
            elif opt_type1 == OptionTypes.AMERICAN_PUT:
                exercise_value = max(k1 - hold_value, 0.0)

            option_vals[index + i_node] = max(exercise_value, hold_value)

    verbose = False

    if verbose:
        print("num_steps1:", num_steps1)
        print("num_steps2:", num_steps2)
        print("u1:", u1, "u2:", u2)
        print("dfs", period_dfs)
        print("probs:", probs)
        print("s:", stock_vals)
        print("v4:", option_vals)

    # We calculate all of the important Greeks in one go
    price = option_vals[0]

    delta = (option_vals[2] - option_vals[1]) / (stock_vals[2] - stock_vals[1])

    delta_up = (option_vals[5] - option_vals[4]) / (stock_vals[5] - stock_vals[4])
    delta_dn = (option_vals[4] - option_vals[3]) / (stock_vals[4] - stock_vals[3])

    gamma = (delta_up - delta_dn) / (stock_vals[2] - stock_vals[1])
    theta = (option_vals[4] - option_vals[0]) / (2.0 * dt1)
    results = np.array([price, delta, gamma, theta])
    return results


########################################################################################



def implied_stock_price(
    s0: float,
    tc: float,
    tu: float,
    kc: float,
    ku: float,
    opt_type_u: OptionTypes,
    r: float,
    q: float,
    vol: float,
) -> float:

    argtuple = (
        opt_type_u,
        tu,
        ku,
        tc,
        r,
        q,
        vol,
        kc,
    )

    sigma = optimize.newton(
        _f,
        x0=s0,
        args=argtuple,
        tol=1e-8,
        maxiter=50,
        fprime2=None,
    )

    return sigma


########################################################################################


def equity_compound_option_bs(
    c_opt_type: OptionTypes,
    u_opt_type: OptionTypes,
    tc: float,
    tu: float,
    kc: float,
    ku: float,
    s0: float,
    ru: float,
    qu: float,
    volatility: float,
    num_steps: int = 200,
) -> float:
    """Value the compound option using an analytical approach if it is
    entirely European style. Otherwise use a Tree approach to handle the
    early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
    Rubinstein (1991). See also Haug page 132."""

    # If the option has any American feature then use the tree
    if (
        c_opt_type == OptionTypes.AMERICAN_CALL
        or u_opt_type == OptionTypes.AMERICAN_CALL
        or c_opt_type == OptionTypes.AMERICAN_PUT
        or u_opt_type == OptionTypes.AMERICAN_PUT
    ):

        v = equity_compound_option_value_tree(
            c_opt_type, u_opt_type, tc, tu, kc, ku, s0, ru, qu, volatility, num_steps
        )

        return v[0]

    # CHECK INTEREST RATES AND IF THERE SHOULD BE TWO RU AND RC ?????
    tc = np.maximum(tc, G_SMALL)
    tu = np.maximum(tc, tu)
    v = np.maximum(volatility, G_SMALL)

    sstar = implied_stock_price(s0, tc, tu, kc, ku, u_opt_type, ru, qu, volatility)

    a1 = (log(s0 / sstar) + (ru - qu + (v**2) / 2.0) * tc) / v / sqrt(tc)
    a2 = a1 - v * sqrt(tc)
    b1 = (log(s0 / ku) + (ru - qu + (v**2) / 2.0) * tu) / v / sqrt(tu)
    b2 = b1 - v * sqrt(tu)

    dqu = exp(-qu * tu)
    dfc = exp(-ru * tc)
    dfu = exp(-ru * tu)
    c = sqrt(tc / tu)

    # Taken from Hull Page 532 (6th edition)

    call_type = OptionTypes.EUROPEAN_CALL
    put_type = OptionTypes.EUROPEAN_PUT

    if c_opt_type == call_type and u_opt_type == call_type:
        v = (
            s0 * dqu * phi2(a1, b1, c)
            - ku * dfu * phi2(a2, b2, c)
            - dfc * kc * normcdf(a2)
        )
    elif c_opt_type == put_type and u_opt_type == call_type:
        v = (
            ku * dfu * phi2(-a2, b2, -c)
            - s0 * dqu * phi2(-a1, b1, -c)
            + dfc * kc * normcdf(-a2)
        )
    elif c_opt_type == call_type and u_opt_type == put_type:
        v = (
            ku * dfu * phi2(-a2, -b2, c)
            - s0 * dqu * phi2(-a1, -b1, c)
            - dfc * kc * normcdf(-a2)
        )
    elif c_opt_type == put_type and u_opt_type == put_type:
        v = (
            s0 * dqu * phi2(a1, -b1, -c)
            - ku * dfu * phi2(a2, -b2, -c)
            + dfc * kc * normcdf(a2)
        )
    else:
        raise FinError("Unknown option type")

    return v


########################################################################################


def equity_compound_option_value_tree(
    c_opt_type: OptionTypes,
    u_opt_type: OptionTypes,
    tc: float,
    tu: float,
    kc: float,
    ku: float,
    s0: float,
    ru: float,
    qu: float,
    volatility: float,
    num_steps: int = 200,
) -> np.ndarray:
    """This function is called if the option has American features."""

    v1 = value_cmpd_once(
        s0,
        ru,
        qu,
        volatility,
        tc,
        tu,
        c_opt_type,
        u_opt_type,
        kc,
        ku,
        num_steps,
    )

    return v1
