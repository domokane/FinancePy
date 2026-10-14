from typing import Any
from math import exp, log, sqrt

import numpy as np
from scipy import optimize
from numba import njit

from ..utils.math import normcdf, phi2
from ..utils.global_vars import G_SMALL
from ..utils.error import FinError
from ..utils.global_types import OptionTypes

from .bs_european_analytic import european_value

########################################################################################


@njit(fastmath=True, cache=True)
def _f(
    stock_price,
    opt_type_u,
    t1,
    t2,
    strike_u,
    r,
    q,
    vol,
    strike_c,
):
    """Root function for the critical stock price.

    The critical stock price S* is defined by

        V_underlying(S*, T1) = K1

    where V_underlying is the value, at the compound-option expiry,
    of the underlying European option.
    """

    if stock_price <= 0.0:
        raise FinError(
            "Stock price must be positive when solving "
            "for the compound-option exercise boundary."
        )

    tau = t2 - t1

    opt_value = european_value(
        stock_price,
        tau,
        strike_u,
        r,
        q,
        vol,
        opt_type_u,
    )

    return opt_value - strike_c

####################################################################################


def implied_stock_price(
    t_c,
    t_u,
    k_c,
    k_u,
    opt_type_u,
    r,
    q,
    vol,
):
    """Calculate the critical stock price for a compound option.

    The critical stock price S* satisfies

        V_underlying(S*, t_u - t_c) = k_c.

    It depends on the compound strike, underlying strike, remaining
    mat_urity, rates, dividend yield, volatility and underlying option
    type. It does not depend on today's stock price.

    A bracketed Brent solver is used because the critical stock price
    must remain positive.
    """

    if t_u <= t_c:
        raise FinError(
            "Underlying option expiry must be after "
            "compound option expiry."
        )

    if k_c <= 0.0:
        raise FinError(
            "Compound option strike must be positive."
        )

    if k_u <= 0.0:
        raise FinError(
            "Underlying option strike must be positive."
        )

    args = (
        opt_type_u,
        t_c,
        t_u,
        k_u,
        r,
        q,
        vol,
        k_c,
    )

    # The stock price must be strictly positive.
    lower = 1.0e-12

    # Start with an economically sensible upper bound and expand it
    # until the root is bracketed.
    upper = max(
        2.0 * k_u,
        2.0 * k_c,
        1.0,
    )

    f_lower = _f(
        lower,
        *args,
    )

    f_upper = _f(
        upper,
        *args,
    )

    max_upper = 1.0e12

    while f_lower * f_upper > 0.0:

        upper *= 2.0

        if upper > max_upper:
            raise FinError(
                "Unable to bracket critical stock price "
                "for compound option."
            )

        f_upper = _f(
            upper,
            *args,
        )

    sstar = optimize.brentq(
        _f,
        lower,
        upper,
        args=args,
        xtol=1.0e-12,
        rtol=1.0e-12,
        maxiter=100,
    )

    return sstar

########################################################################################


@njit(fastmath=True, cache=True)
def value_cmpd_once(
    s: float,
    r: float,
    q: float,
    volatility: float,
    t1: float,
    t2: float,
    opt_type1: int,
    opt_type2: int,
    k1: float,
    k2: float,
    num_steps: int,
) -> np.ndarray:

    if num_steps < 3:
        num_steps = 3

    # TODO EUROPEAN call-put works but AMERICAN call-put needs to be tested

    # Need equally spaced time intervals for a recombining tree
    # Downside is that we may not measure periods exactly
    dt = t2 / num_steps

    # Map the compound-option expiry to the nearest tree time.
    # Using int() would systematically place the compound expiry
    # before its act_ual date and introduce a pricing bias.
    num_steps1 = int(round(t1 / dt))
    num_steps1 = max(1, min(num_steps1, num_steps - 1))
    num_steps2 = num_steps - num_steps1

    dt1 = dt
    dt2 = dt

    # the number of nodes on the tree
    num_nodes = (num_steps + 1) * (num_steps + 2) // 2
    option_vals = np.zeros(num_nodes)

    u1 = np.exp(volatility * np.sqrt(dt))
    d1 = 1.0 / u1
    u2 = np.exp(volatility * np.sqrt(dt))
    d2 = 1.0 / u2

    p = (np.exp((r - q) * dt) - d1) / (u1 - d1)

    if not np.isfinite(p) or p < 0.0 or p > 1.0:
        raise FinError(
            "Invalid CRR probability; increase num_steps"
        )

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

    for i_node in range(0, num_steps + 1):
        s = stock_vals[index + i_node]
        if (
            opt_type2 == OptionTypes.EUROPEAN_CALL.value
            or opt_type2 == OptionTypes.AMERICAN_CALL.value
        ):
            option_vals[index + i_node] = max(s - k2, 0.0)
        elif (
            opt_type2 == OptionTypes.EUROPEAN_PUT.value
            or opt_type2 == OptionTypes.AMERICAN_PUT.value
        ):
            option_vals[index + i_node] = max(k2 - s, 0.0)

    # begin backward steps from expiry at t2 to first expiry at time t1
    for i_time in range(num_steps - 1, num_steps1 - 1, -1):
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

            if opt_type2 == OptionTypes.AMERICAN_CALL.value:
                exercise_value = max(s - k2, 0.0)
            elif opt_type2 == OptionTypes.AMERICAN_PUT.value:
                exercise_value = max(k2 - s, 0.0)

            option_vals[index + i_node] = max(exercise_value, hold_value)

    # Now do payoff at the end of the first expiry period at t1
    i_time = num_steps1
    index = i_time * (i_time + 1) // 2

    for i_node in range(0, i_time + 1):
        #     s = stock_vals[index + i_node]
        #     next_index = (i_time + 1) * (i_time + 2) // 2
        #     next_node_dn = next_index + i_node
        #     next_node_up = next_index + i_node + 1
        #     v_up = option_vals[next_node_up]
        #     v_dn = option_vals[next_node_dn]
        #     future_exp_val = probs[i_time] * v_up
        #     future_exp_val += (1.0 - probs[i_time]) * v_dn
        #     hold_value = period_dfs[i_time] * fut_ure_exp_val

        underlying_value = option_vals[index + i_node]

        if (
            opt_type1 == OptionTypes.EUROPEAN_CALL.value
            or opt_type1 == OptionTypes.AMERICAN_CALL.value
        ):
            option_vals[index + i_node] = max(underlying_value - k1, 0.0)
        elif (
            opt_type1 == OptionTypes.EUROPEAN_PUT.value
            or opt_type1 == OptionTypes.AMERICAN_PUT.value
        ):
            option_vals[index + i_node] = max(k1 - underlying_value, 0.0)

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

            if opt_type1 == OptionTypes.AMERICAN_CALL.value:
                exercise_value = max(hold_value - k1, 0.0)
            elif opt_type1 == OptionTypes.AMERICAN_PUT.value:
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

    gamma = 2.0 * (delta_up - delta_dn) / (stock_vals[5] - stock_vals[3])
    theta = (option_vals[4] - option_vals[0]) / (2.0 * dt1)
    results = np.array([price, delta, gamma, theta])
    return results


########################################################################################

def equity_compound_option_bs(
    c_opt_type: int,
    u_opt_type: int,
    t_c: float,
    t_u: float,
    k_c: float,
    k_u: float,
    s0: float,
    r_c: float,
    r_u: float,
    q_c: float,
    q_u: float,
    volatility: float
) -> float:
    """Value the compound option using an analytical approach if it is
    entirely European style. Otherwise use a Tree approach to handle the
    early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
    Rubinstein (1991). See also Haug page 132."""

    call_type = OptionTypes.EUROPEAN_CALL.value
    put_type = OptionTypes.EUROPEAN_PUT.value

    if not np.isfinite(t_c) or not np.isfinite(t_u):
        raise FinError("Maturities must be finite")

    if t_c <= 0.0 or t_u <= t_c:
        raise FinError(
            "Require 0 < compound expiry < underlying expiry"
        )

    if not np.isfinite(k_c) or k_c < 0.0:
        raise FinError("Invalid compound strike")

    if not np.isfinite(k_u) or k_u <= 0.0:
        raise FinError("Invalid underlying strike")

    if not np.isfinite(s0) or s0 <= 0.0:
        raise FinError("Invalid stock price")

    if not np.isfinite(r_u) or not np.isfinite(q_u):
        raise FinError("Invalid underlying option interest rate or dividend yield")

    if not np.isfinite(r_c) or not np.isfinite(q_c):
        raise FinError("Invalid compoubnd option interest rate or dividend yield")

    if not np.isfinite(volatility) or volatility < 0.0:
        raise FinError("Invalid volatility")

    v = max(volatility, G_SMALL)

    tau = t_u - t_c

    r_cu = (r_u * t_u - r_c * t_c) / tau
    q_cu = (q_u * t_u - q_c * t_c) / tau

    underlying_value = european_value(
        s0, t_u, k_u, r_u, q_u, v, u_opt_type
    )

    if k_c == 0.0:
        if c_opt_type == call_type:
            return underlying_value
        if c_opt_type == put_type:
            return 0.0

    if u_opt_type == put_type:
        max_put_value = k_u * exp(-r_cu * (t_u - t_c))

        if k_c >= max_put_value:
            if c_opt_type == call_type:
                return 0.0
            if c_opt_type == put_type:
                return k_c * exp(-r_c * t_c) - underlying_value

    sstar = implied_stock_price(
        t_c, t_u, k_c, k_u,
        u_opt_type, r_cu, q_cu, v
    )

    a1 = (log(s0 / sstar) + (r_c - q_c + 0.50 * (v**2)) * t_c) / v / sqrt(t_c)
    a2 = a1 - v * sqrt(t_c)
    b1 = (log(s0 / k_u) + (r_u - q_u + 0.50 * (v**2)) * t_u) / v / sqrt(t_u)
    b2 = b1 - v * sqrt(t_u)

    dfc = exp(-r_c * t_c)
    dfu = exp(-r_u * t_u)

    dqu = exp(-q_u * t_u)
    c = sqrt(t_c / t_u)

    # Taken from Hull Page 532 (6th edition)

    call_type = OptionTypes.EUROPEAN_CALL.value
    put_type = OptionTypes.EUROPEAN_PUT.value

    if c_opt_type == call_type and u_opt_type == call_type:
        v = (
            s0 * dqu * phi2(a1, b1, c)
            - k_u * dfu * phi2(a2, b2, c)
            - dfc * k_c * normcdf(a2)
        )
    elif c_opt_type == put_type and u_opt_type == call_type:
        v = (
            k_u * dfu * phi2(-a2, b2, -c)
            - s0 * dqu * phi2(-a1, b1, -c)
            + dfc * k_c * normcdf(-a2)
        )
    elif c_opt_type == call_type and u_opt_type == put_type:
        v = (
            k_u * dfu * phi2(-a2, -b2, c)
            - s0 * dqu * phi2(-a1, -b1, c)
            - dfc * k_c * normcdf(-a2)
        )
    elif c_opt_type == put_type and u_opt_type == put_type:
        v = (
            s0 * dqu * phi2(a1, -b1, -c)
            - k_u * dfu * phi2(a2, -b2, -c)
            + dfc * k_c * normcdf(a2)
        )
    else:
        raise FinError("Unknown option type")

    return v


########################################################################################

@njit(fastmath=True, cache=True)
def equity_compound_option_value_tree(
    c_opt_type: int,
    u_opt_type: int,
    t_c: float,
    t_u: float,
    k_c: float,
    k_u: float,
    s0: float,
    ru: float,
    qu: float,
    volatility: float,
    num_steps: int = 200,
) -> np.ndarray:
    """This function is called if the option has American feat_ures."""

    v1 = value_cmpd_once(
        s0,
        ru,
        qu,
        volatility,
        t_c,
        t_u,
        c_opt_type,
        u_opt_type,
        k_c,
        k_u,
        num_steps,
    )

    v2 = value_cmpd_once(
        s0,
        ru,
        qu,
        volatility,
        t_c,
        t_u,
        c_opt_type,
        u_opt_type,
        k_c,
        k_u,
        num_steps+1,
    )

    v = 0.5*(v1+v2)

    return v
