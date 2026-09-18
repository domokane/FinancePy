##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO - MUST ADD ACCRUED INTEREST TO MODEL!!!!

from math import exp, sqrt

from numba import njit
import numpy as np

from ..utils.error import FinError
from ..market.curves.interpolator import InterpTypes, _uinterpolate

########################################################################################


@njit(fastmath=True, cache=True)
def value_convertible(
    t_mat,
    face_amount,
    cpn_times,
    cpn_flows,
    call_times,
    call_prices,
    put_times,
    put_prices,
    conv_ratio,
    start_convert_time,
    # Market inputs
    stock_price,
    df_times,
    df_values,
    dividend_times,
    dividend_yields,
    stock_volatility,
    credit_spread,
    recovery_rate,
    # Tree details
    num_steps_per_year,
):

    interp = InterpTypes.FLAT_FWD_RATES.value

    if len(cpn_times) > 0:
        if cpn_times[-1] > t_mat:
            raise FinError("Coupon after maturity")

    if len(call_times) > 0:
        if call_times[-1] > t_mat:
            raise FinError("Call times after maturity")

    if len(put_times) > 0:
        if put_times[-1] > t_mat:
            raise FinError("Put times after maturity")

    if len(df_times) > 0:
        if df_times[-1] > t_mat:
            raise FinError("Discount times after maturity")

    if len(dividend_times) > 0:
        if dividend_times[-1] > t_mat:
            raise FinError("Dividend times after maturity")

    if credit_spread < 0.0:
        raise FinError("Credit spread negative.")

    if recovery_rate < 0.0 or recovery_rate > 1.0:
        raise FinError("Recovery rate should be between 0 and 1.")

    if stock_volatility < 0.0:
        raise FinError("Stock volatility cannot be negative.")

    if num_steps_per_year < 1:
        raise FinError("Num Steps per year must more than 1.")

    if len(dividend_times) > 0.0:
        if dividend_times[-1] > t_mat:
            raise FinError("Last dividend is after bond maturity.")

    if recovery_rate > 0.999 or recovery_rate < 0.0:
        raise FinError("Recovery rate must be between 0 and 0.999.")

    num_times = int(num_steps_per_year * t_mat) + 1  # add one for today time 0

    if num_times < 5:
        raise FinError("Numsteps must be greater than 5.")

    num_levels = num_times

    # this is the size of the step
    dt = t_mat / (num_times - 1)

    tree_times = np.linspace(0.0, t_mat, num_times)
    tree_dfs = np.zeros(num_times)
    for i in range(0, num_times):
        df = _uinterpolate(tree_times[i], df_times, df_values, interp)
        tree_dfs[i] = df

    h = credit_spread / (1.0 - recovery_rate)
    survival_prob = exp(-h * dt)

    # map coupons onto tree but preserve their present value using risky dfs
    tree_flows = np.zeros(num_times)
    num_cpns = len(cpn_times)
    for i in range(0, num_cpns):
        flow_time = cpn_times[i]
        n = int(round(flow_time / dt, 0))
        tree_time = tree_times[n]
        df_flow = _uinterpolate(flow_time, df_times, df_values, interp)
        df_flow *= exp(-h * flow_time)
        df_tree = _uinterpolate(tree_time, df_times, df_values, interp)
        df_tree *= exp(-h * tree_time)
        tree_flows[n] += cpn_flows[i] * 1.0 * df_flow / df_tree

    # map call onto tree - must have no calls at high value
    tree_call_value = np.ones(num_times) * face_amount * 1000.0
    num_calls = len(call_times)
    for i in range(0, num_calls):
        call_time = call_times[i]
        n = int(round(call_time / dt, 0))
        tree_call_value[n] = call_prices[i]

    # map puts onto tree
    tree_put_value = np.zeros(num_times)
    num_puts = len(put_times)
    for i in range(0, num_puts):
        put_time = put_times[i]
        n = int(round(put_time / dt, 0))
        tree_put_value[n] = put_prices[i]

    # map discrete dividend yields onto tree dates when they are made
    tree_dividend_yld = np.zeros(num_times)
    num_dividends = len(dividend_times)
    for i in range(0, num_dividends):
        dividend_time = dividend_times[i]
        n = int(round(dividend_time / dt, 0))
        tree_dividend_yld[n] = dividend_yields[i]

    # Set up the tree of stock prices using a 2D matrix - half the matrix is
    # unused but this may be a cost worth bearing for simpler code. Review.
    tree_stock_value = np.zeros(shape=(num_times, num_levels))
    e = stock_volatility**2 - h
    if e < 0.0:
        raise FinError("Volatility squared minus the hazard rate is negative.")

    u = exp(sqrt(e * dt))
    d = 1.0 / u
    u2 = u * u
    tree_stock_value[0, 0] = stock_price
    for i_time in range(1, num_times):
        s = tree_stock_value[i_time - 1, 0] * d
        tree_stock_value[i_time, 0] = s

        for i_node in range(1, i_time + 1):
            s = s * u2
            tree_stock_value[i_time, i_node] = s

        # we now reduce all stocks by the same yield amount at the same date
        y = tree_dividend_yld[i_time]
        for i_node in range(0, i_time + 1):
            tree_stock_value[i_time, i_node] *= 1.0 - y

    # set up the tree of conversion values. Before allowed to convert the
    # conversion value must be set equal to zero

    tree_convert_value = np.zeros(shape=(num_times, num_levels))
    for i_time in range(0, num_times):
        if tree_times[i_time] >= start_convert_time:
            for i_node in range(0, i_time + 1):
                s = tree_stock_value[i_time, i_node]
                tree_convert_value[i_time, i_node] = s * conv_ratio * 1.0

    #    print_tree(tree_convert_value)

    tree_convert_bond_value = np.zeros(shape=(num_times, num_levels))

    # store probability of up move as a function of time on the tree
    tree_probs_up = np.zeros(num_times)
    tree_probs_dn = np.zeros(num_times)
    q = 0.0  # we have discrete dividends paid as dividend yields only
    for i_time in range(1, num_times):
        a = tree_dfs[i_time - 1] / tree_dfs[i_time] * exp(-q * dt)
        tree_probs_up[i_time] = (a - d * survival_prob) / (u - d)
        tree_probs_dn[i_time] = (u * survival_prob - a) / (u - d)
    #        r = log(a)/dt
    #        n_min = r*r / stock_volatility / stock_volatility

    if np.any(tree_probs_up > 1.0):
        raise FinError("p_up > 1.0. Increase time steps.")

    ###########################################################################
    # work backwards by first setting values at bond maturity date
    ###########################################################################

    flow = tree_flows[num_times - 1]
    bullet_pv = (1.0 + flow) * face_amount
    for i_node in range(0, num_levels):
        conv_value = tree_convert_value[num_times - 1, i_node]
        tree_convert_bond_value[num_times - 1, i_node] = max(bullet_pv, conv_value)

    #  begin backward steps from expiry
    for i_time in range(num_times - 2, -1, -1):

        p_up = tree_probs_up[i_time + 1]
        p_dn = tree_probs_dn[i_time + 1]
        p_def = 1.0 - survival_prob
        df = tree_dfs[i_time + 1] / tree_dfs[i_time]
        call = tree_call_value[i_time]
        put = tree_put_value[i_time]
        flow = tree_flows[i_time]

        for i_node in range(0, i_time + 1):
            fut_value_up = tree_convert_bond_value[i_time + 1, i_node + 1]
            fut_value_dn = tree_convert_bond_value[i_time + 1, i_node]
            hold = p_up * fut_value_up + p_dn * fut_value_dn  # p_up already embeds Q
            hold_pv = df * hold + p_def * df * recovery_rate * face_amount + flow * face_amount
            conv = tree_convert_value[i_time, i_node]
            value = min(max(hold_pv, conv, put), call)
            tree_convert_bond_value[i_time, i_node] = value

        bullet_pv = df * bullet_pv * survival_prob
        bullet_pv += p_def * df * recovery_rate * face_amount
        bullet_pv += flow * face_amount

    price = tree_convert_bond_value[0, 0]
    delta = (tree_convert_bond_value[1, 1] - tree_convert_bond_value[1, 0]) / (
        tree_stock_value[1, 1] - tree_stock_value[1, 0]
    )
    delta_up = (tree_convert_bond_value[2, 3] - tree_convert_bond_value[2, 2]) / (
        tree_stock_value[2, 3] - tree_stock_value[2, 2]
    )
    delta_dn = (tree_convert_bond_value[2, 2] - tree_convert_bond_value[2, 1]) / (
        tree_stock_value[2, 2] - tree_stock_value[2, 1]
    )
    gamma = (delta_up - delta_dn) / (tree_stock_value[1, 1] - tree_stock_value[1, 0])
    theta = (tree_convert_bond_value[2, 2] - tree_convert_bond_value[0, 0]) / (2.0 * dt)
    results = np.array([price, bullet_pv, delta, gamma, theta])
    return results
