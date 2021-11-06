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
                 ccInterestRate,  # continuously compounded
                 ccDividendRate,  # continuously compounded
                 volatility,  # Black scholes volatility
                 num_steps_per_year,
                 time_to_expiry,
                 option_type,
                 strike_price,
                 isEven):
    """ Value an American option using a Binomial Treee """

    num_steps = int(num_steps_per_year * time_to_expiry)

    if num_steps < 30:
        num_steps = 30

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
    r = ccInterestRate
    q = ccDividendRate

    # the number of nodes on the tree
    num_nodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stock_values = np.zeros(num_nodes)
    stock_values[0] = stock_price

    option_values = np.zeros(num_nodes)
    u = np.exp(volatility * np.sqrt(dt))
    d = 1.0 / u
    sLow = stock_price

    probs = np.zeros(num_steps)
    periodDiscountFactors = np.zeros(num_steps)

    # store time independent information for later use in tree
    for iTime in range(0, num_steps):
        a = np.exp((r - q) * dt)
        probs[iTime] = (a - d) / (u - d)
        periodDiscountFactors[iTime] = np.exp(-r * dt)

    for iTime in range(1, num_steps + 1):
        sLow *= d
        s = sLow
        for iNode in range(0, iTime + 1):
            index = 0.5 * iTime * (iTime + 1)
            stock_values[int(index + iNode)] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * num_steps * (num_steps + 1))

    for iNode in range(0, iTime + 1):

        s = stock_values[index + iNode]

        if option_type == OptionTypes.EUROPEAN_CALL.value:
            option_values[index + iNode] = np.maximum(s - strike_price, 0.0)
        elif option_type == OptionTypes.EUROPEAN_PUT.value:
            option_values[index + iNode] = np.maximum(strike_price - s, 0.0)
        elif option_type == OptionTypes.AMERICAN_CALL.value:
            option_values[index + iNode] = np.maximum(s - strike_price, 0.0)
        elif option_type == OptionTypes.AMERICAN_PUT.value:
            option_values[index + iNode] = np.maximum(strike_price - s, 0.0)

    # begin backward steps from expiry to value date
    for iTime in range(num_steps - 1, -1, -1):

        index = int(0.5 * iTime * (iTime + 1))

        for iNode in range(0, iTime + 1):

            s = stock_values[index + iNode]

            exerciseValue = 0.0

            if option_type == OptionTypes.EUROPEAN_CALL.value:
                exerciseValue = 0.0
            elif option_type == OptionTypes.EUROPEAN_PUT.value:
                exerciseValue = 0.0
            elif option_type == OptionTypes.AMERICAN_CALL.value:
                exerciseValue = np.maximum(s - strike_price, 0.0)
            elif option_type == OptionTypes.AMERICAN_PUT.value:
                exerciseValue = np.maximum(strike_price - s, 0.0)

            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))

            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1

            vUp = option_values[nextNodeUp]
            vDn = option_values[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            if option_type == OptionTypes.EUROPEAN_CALL.value:
                option_values[index + iNode] = holdValue
            elif option_type == OptionTypes.EUROPEAN_PUT.value:
                option_values[index + iNode] = holdValue
            elif option_type == OptionTypes.AMERICAN_CALL.value:
                option_values[index +
                              iNode] = np.maximum(exerciseValue, holdValue)
            elif option_type == OptionTypes.AMERICAN_PUT.value:
                option_values[index +
                              iNode] = np.maximum(exerciseValue, holdValue)

    # We calculate all of the important Greeks in one go
    price = option_values[0]
    delta = (option_values[2] - option_values[1]) / \
        (stock_values[2] - stock_values[1])
    deltaUp = (option_values[5] - option_values[4]) / \
        (stock_values[5] - stock_values[4])
    deltaDn = (option_values[4] - option_values[3]) / \
        (stock_values[4] - stock_values[3])
    gamma = (deltaUp - deltaDn) / (stock_values[2] - stock_values[1])
    theta = (option_values[4] - option_values[0]) / (2.0 * dt)
    results = np.array([price, delta, gamma, theta])
    return results

###############################################################################


def crr_tree_val_avg(stock_price,
                     ccInterestRate,  # continuously compounded
                     ccDividendRate,  # continuously compounded
                     volatility,  # Black scholes volatility
                     num_steps_per_year,
                     time_to_expiry,
                     option_type,
                     strike_price):
    """ Calculate the average values off the tree using an even and an odd
    number of time steps. """

    value1 = crr_tree_val(stock_price,
                          ccInterestRate,
                          ccDividendRate,
                          volatility,
                          num_steps_per_year,
                          time_to_expiry,
                          option_type,
                          strike_price,
                          1)  # even

    value2 = crr_tree_val(stock_price,
                          ccInterestRate,
                          ccDividendRate,
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
