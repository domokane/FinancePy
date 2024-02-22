##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ..utils.global_types import OptionTypes

import numpy as np
from numba import njit, float64, int64

bump = 1e-4

###############################################################################
###############################################################################


@njit(float64[:](float64, float64, float64, float64, int64, float64, int64,
                 float64, int64), fastmath=True, cache=True)
def crr_tree_val(stock_price,
                 interest_rate,  # continuously compounded
                 dividend_rate,  # continuously compounded
                 volatility,  # Black scholes volatility
                 num_steps_per_year,
                 time_to_expiry,
                 option_type,
                 strike_price,
                 isEven):
    """ Value an American option using a Binomial Treee """

    num_steps = int(num_steps_per_year * time_to_expiry)
    num_steps = max(num_steps, 30)

    # OVERRIDE JUST TO SEE
    num_steps = num_steps_per_year

    # if the number of steps is even but we want odd then make it odd
    if num_steps % 2 == 0 and isEven == 0:
        num_steps += 1
    elif num_steps % 2 == 1 and isEven == 1:
        num_steps += 1

#    print(num_steps)
    # this is the size of the step
    dt = time_to_expiry / num_steps
    r = interest_rate
    q = dividend_rate

    # the number of nodes on the tree
    num_nodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stock_values = np.zeros(num_nodes)
    stock_values[0] = stock_price

    option_values = np.zeros(num_nodes)
    u = np.exp(volatility * np.sqrt(dt))
    d = 1.0 / u
    s_low = stock_price

    probs = np.zeros(num_steps)
    period_dfs = np.zeros(num_steps)

    # store time independent information for later use in tree
    for i_time in range(0, num_steps):
        a = np.exp((r - q) * dt)
        probs[i_time] = (a - d) / (u - d)
        period_dfs[i_time] = np.exp(-r * dt)

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

        if option_type == OptionTypes.EUROPEAN_CALL.value:
            option_values[index + i_node] = np.maximum(s - strike_price, 0.0)
        elif option_type == OptionTypes.EUROPEAN_PUT.value:
            option_values[index + i_node] = np.maximum(strike_price - s, 0.0)
        elif option_type == OptionTypes.AMERICAN_CALL.value:
            option_values[index + i_node] = np.maximum(s - strike_price, 0.0)
        elif option_type == OptionTypes.AMERICAN_PUT.value:
            option_values[index + i_node] = np.maximum(strike_price - s, 0.0)

    # begin backward steps from expiry to value date
    for i_time in range(num_steps - 1, -1, -1):

        index = int(0.5 * i_time * (i_time + 1))

        for i_node in range(0, i_time + 1):

            s = stock_values[index + i_node]

            exercise_value = 0.0

            if option_type == OptionTypes.EUROPEAN_CALL.value:
                exercise_value = 0.0
            elif option_type == OptionTypes.EUROPEAN_PUT.value:
                exercise_value = 0.0
            elif option_type == OptionTypes.AMERICAN_CALL.value:
                exercise_value = np.maximum(s - strike_price, 0.0)
            elif option_type == OptionTypes.AMERICAN_PUT.value:
                exercise_value = np.maximum(strike_price - s, 0.0)

            next_index = int(0.5 * (i_time + 1) * (i_time + 2))

            next_node_dn = next_index + i_node
            next_node_up = next_index + i_node + 1

            v_up = option_values[next_node_up]
            v_dn = option_values[next_node_dn]
            future_expected_value = probs[i_time] * v_up
            future_expected_value += (1.0 - probs[i_time]) * v_dn
            hold_value = period_dfs[i_time] * future_expected_value

            if option_type == OptionTypes.EUROPEAN_CALL.value:
                option_values[index + i_node] = hold_value
            elif option_type == OptionTypes.EUROPEAN_PUT.value:
                option_values[index + i_node] = hold_value
            elif option_type == OptionTypes.AMERICAN_CALL.value:
                option_values[index +
                              i_node] = np.maximum(exercise_value, hold_value)
            elif option_type == OptionTypes.AMERICAN_PUT.value:
                option_values[index +
                              i_node] = np.maximum(exercise_value, hold_value)

    # We calculate all of the important Greeks in one go
    price = option_values[0]
    delta = (option_values[2] - option_values[1]) / \
        (stock_values[2] - stock_values[1])
    delta_up = (option_values[5] - option_values[4]) / \
        (stock_values[5] - stock_values[4])
    delta_dn = (option_values[4] - option_values[3]) / \
        (stock_values[4] - stock_values[3])
    gamma = (delta_up - delta_dn) / (stock_values[2] - stock_values[1])
    theta = (option_values[4] - option_values[0]) / (2.0 * dt)
    results = np.array([price, delta, gamma, theta])
    return results

###############################################################################


def crr_tree_val_avg(stock_price,
                     interest_rate,  # continuously compounded
                     dividend_rate,  # continuously compounded
                     volatility,  # Black scholes volatility
                     num_steps_per_year,
                     time_to_expiry,
                     option_type,
                     strike_price):
    """ Calculate the average values off the tree using an even and an odd
    number of time steps. """

    value1 = crr_tree_val(stock_price,
                          interest_rate,
                          dividend_rate,
                          volatility,
                          num_steps_per_year,
                          time_to_expiry,
                          option_type,
                          strike_price,
                          1)  # even

    value2 = crr_tree_val(stock_price,
                          interest_rate,
                          dividend_rate,
                          volatility,
                          num_steps_per_year,
                          time_to_expiry,
                          option_type,
                          strike_price,
                          0)  # odd

    v = (value1 + value2) / 2.0
    res = {'value': v[0], 'delta': v[1], 'gamma': v[2], 'theta': v[3]}
    return res

###############################################################################
