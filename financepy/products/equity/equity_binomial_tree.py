##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
from enum import Enum

import numpy as np
from numba import njit, float64, int64

from ...utils.error import FinError
from ...utils.global_vars import g_days_in_year
from ...utils.math import heaviside

###############################################################################


class EquityTreePayoffTypes(Enum):
    FWD_CONTRACT = 1
    VANILLA_OPTION = 2
    DIGITAL_OPTION = 3
    POWER_CONTRACT = 4
    POWER_OPTION = 5
    LOG_CONTRACT = 6
    LOG_OPTION = 7


class EquityTreeExerciseTypes(Enum):
    EUROPEAN = 1
    AMERICAN = 2


###############################################################################


@njit
def _validate_payoff(payoff_type, payoff_params):

    num_params = 0

    if payoff_type == EquityTreePayoffTypes.FWD_CONTRACT.value:
        num_params = 1
    elif payoff_type == EquityTreePayoffTypes.VANILLA_OPTION.value:
        num_params = 2
    elif payoff_type == EquityTreePayoffTypes.DIGITAL_OPTION.value:
        num_params = 2
    elif payoff_type == EquityTreePayoffTypes.POWER_CONTRACT.value:
        num_params = 2
    elif payoff_type == EquityTreePayoffTypes.POWER_OPTION.value:
        num_params = 3
    elif payoff_type == EquityTreePayoffTypes.LOG_CONTRACT.value:
        num_params = 0
    elif payoff_type == EquityTreePayoffTypes.LOG_OPTION.value:
        num_params = 1
    else:
        raise FinError("Unknown payoff type")

    if len(payoff_params) != num_params:
        raise FinError(
            "Number of parameters required for "
            + str(payoff_type)
            + " must be "
            + str(num_params)
        )

    return None


###############################################################################


@njit(float64(float64, int64, float64[:]), fastmath=True, cache=True)
def _payoff_value(s, payoff_type, payoff_params):

    if payoff_type == EquityTreePayoffTypes.FWD_CONTRACT.value:
        payoff = payoff_params[0] * s
    elif payoff_type == EquityTreePayoffTypes.VANILLA_OPTION.value:
        payoff = max(payoff_params[0] * (s - payoff_params[1]), 0.0)
    elif payoff_type == EquityTreePayoffTypes.DIGITAL_OPTION.value:
        payoff = heaviside(payoff_params[0] * (s - payoff_params[1]))
    elif payoff_type == EquityTreePayoffTypes.POWER_CONTRACT.value:
        payoff = payoff_params[0] * (s ** payoff_params[1])
    elif payoff_type == EquityTreePayoffTypes.POWER_OPTION.value:
        payoff = max(
            payoff_params[0] * ((s ** payoff_params[2]) - payoff_params[1]),
            0.0,
        )
    elif payoff_type == EquityTreePayoffTypes.LOG_CONTRACT.value:
        payoff = log(s)
    elif payoff_type == EquityTreePayoffTypes.LOG_OPTION.value:
        payoff = max(log(s) - payoff_params[0], 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff


###############################################################################


@njit(fastmath=True, cache=True)
def _value_once(
    stock_price,
    r,
    q,
    volatility,
    num_steps,
    time_to_expiry,
    payoff_type,
    exercise_type,
    payoff_params,
):

    num_steps = max(num_steps, 3)

    #        validate_payoff(payoff_type.value,payoff_params)

    payoff_typeValue = payoff_type.value

    # this is the size of the step
    dt = time_to_expiry / num_steps

    # the number of nodes on the tree
    num_nodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stock_values = np.zeros(num_nodes)
    stock_values[0] = stock_price

    option_values = np.zeros(num_nodes)
    u = exp(volatility * sqrt(dt))
    d = 1.0 / u
    s_low = stock_price

    probs = np.zeros(num_steps)
    period_dfs = np.zeros(num_steps)

    # store time independent information for later use in tree
    for i_time in range(0, num_steps):
        a = exp((r - q) * dt)
        probs[i_time] = (a - d) / (u - d)
        period_dfs[i_time] = exp(-r * dt)

    for i_time in range(1, num_steps + 1):
        s_low *= d
        s = s_low
        for i_node in range(0, i_time + 1):
            index = 0.5 * i_time * (i_time + 1)
            stock_values[int(index + i_node)] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * num_steps * (num_steps + 1))

    for i_node in range(0, i_time + 1):
        s = stock_values[index + i_node]
        option_values[index + i_node] = _payoff_value(
            s, payoff_typeValue, payoff_params
        )

    # begin backward steps from expiry
    for i_time in range(num_steps - 1, -1, -1):

        index = int(0.5 * i_time * (i_time + 1))

        for i_node in range(0, i_time + 1):

            next_index = int(0.5 * (i_time + 1) * (i_time + 2))
            next_node_dn = next_index + i_node
            next_node_up = next_node_dn + 1
            v_up = option_values[next_node_up]
            v_dn = option_values[next_node_dn]
            future_expected_value = probs[i_time] * v_up
            future_expected_value += (1.0 - probs[i_time]) * v_dn
            hold_value = period_dfs[i_time] * future_expected_value

            if exercise_type == EquityTreeExerciseTypes.EUROPEAN:
                option_values[index + i_node] = hold_value
            elif exercise_type == EquityTreeExerciseTypes.AMERICAN:
                s = stock_values[index + i_node]
                exercise_value = _payoff_value(
                    s, payoff_typeValue, payoff_params
                )
                option_values[index + i_node] = max(exercise_value, hold_value)

    price = option_values[0]
    delta = (option_values[2] - option_values[1]) / (
        stock_values[2] - stock_values[1]
    )
    delta_up = (option_values[5] - option_values[4]) / (
        stock_values[5] - stock_values[4]
    )
    delta_dn = (option_values[4] - option_values[3]) / (
        stock_values[4] - stock_values[3]
    )
    gamma = (delta_up - delta_dn) / (stock_values[2] - stock_values[1])
    theta = (option_values[4] - option_values[0]) / (2.0 * dt)
    results = np.array([price, delta, gamma, theta])
    return results


###############################################################################


class EquityBinomialTree:

    def __init__(self):
        pass

    #        self.m_option_values = np.zeros()
    #        self.m_stock_values = np.zeros()
    #        self.m_upProbabilities = np.zeros()
    #
    #       self.m_num_steps = 10
    #        self.m_num_nodes = 10

    ###############################################################################

    def value(
        self,
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        value_dt,
        payoff,
        expiry_dt,
        payoff_type,
        exercise_type,
        payoff_params,
    ):

        # do some validation
        t_exp = (expiry_dt - value_dt) / g_days_in_year
        r = discount_curve.zero_rate(expiry_dt)

        dq = dividend_curve.df(expiry_dt)
        q = -np.log(dq) / t_exp

        price1 = _value_once(
            stock_price,
            r,
            q,
            volatility,
            num_steps,
            t_exp,
            payoff_type,
            exercise_type,
            payoff_params,
        )

        # Can I reuse the same tree ?
        price2 = _value_once(
            stock_price,
            r,
            q,
            volatility,
            num_steps + 1,
            t_exp,
            payoff_type,
            exercise_type,
            payoff_params,
        )

        price = (price1 + price2) / 2.0

        return price


###############################################################################
