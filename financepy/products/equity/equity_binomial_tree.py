##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
from enum import Enum

import numpy as np
from numba import jit, njit, float64, int64

from ...utils.error import FinError
from ...utils.global_vars import gDaysInYear
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


@jit
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
            "Number of parameters required for " +
            str(payoff_type) +
            " must be " +
            str(num_params))

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
        payoff = payoff_params[0] * (s**payoff_params[1])
    elif payoff_type == EquityTreePayoffTypes.POWER_OPTION.value:
        payoff = max(payoff_params[0] *
                     ((s**payoff_params[2]) -
                      payoff_params[1]), 0.0)
    elif payoff_type == EquityTreePayoffTypes.LOG_CONTRACT.value:
        payoff = log(s)
    elif payoff_type == EquityTreePayoffTypes.LOG_OPTION.value:
        payoff = max(log(s) - payoff_params[0], 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff

###############################################################################


@njit(fastmath=True, cache=True)
def _value_once(stock_price,
                r,
                q,
                volatility,
                num_steps,
                time_to_expiry,
                payoff_type,
                exercise_type,
                payoff_params):

    if num_steps < 3:
        num_steps = 3

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
    sLow = stock_price

    probs = np.zeros(num_steps)
    periodDiscountFactors = np.zeros(num_steps)

    # store time independent information for later use in tree
    for iTime in range(0, num_steps):
        a = exp((r - q) * dt)
        probs[iTime] = (a - d) / (u - d)
        periodDiscountFactors[iTime] = exp(-r * dt)

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
        option_values[index +
                      iNode] = _payoff_value(s, payoff_typeValue, payoff_params)

    # begin backward steps from expiry
    for iTime in range(num_steps - 1, -1, -1):

        index = int(0.5 * iTime * (iTime + 1))

        for iNode in range(0, iTime + 1):

            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextNodeDn + 1
            vUp = option_values[nextNodeUp]
            vDn = option_values[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            if exercise_type == EquityTreeExerciseTypes.EUROPEAN:
                option_values[index + iNode] = holdValue
            elif exercise_type == EquityTreeExerciseTypes.AMERICAN:
                s = stock_values[index + iNode]
                exerciseValue = _payoff_value(
                    s, payoff_typeValue, payoff_params)
                option_values[index + iNode] = max(exerciseValue, holdValue)

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


class EquityBinomialTree():

    def __init__(self):
        pass
#        self.m_option_values = np.zeros()
#        self.m_stock_values = np.zeros()
#        self.m_upProbabilities = np.zeros()
#
#       self.m_num_steps = 10
#        self.m_num_nodes = 10

###############################################################################

    def value(self,
              stock_price,
              discount_curve,
              dividend_curve,
              volatility,
              num_steps,
              valuation_date,
              payoff,
              expiry_date,
              payoff_type,
              exercise_type,
              payoff_params):

        # do some validation
        texp = (expiry_date - valuation_date) / gDaysInYear
        r = discount_curve.zero_rate(expiry_date)

        dq = dividend_curve.df(expiry_date)
        q = -np.log(dq)/texp

        price1 = _value_once(stock_price,
                             r,
                             q,
                             volatility,
                             num_steps,
                             texp,
                             payoff_type,
                             exercise_type,
                             payoff_params)

        # Can I reuse the same tree ?
        price2 = _value_once(stock_price,
                             r,
                             q,
                             volatility,
                             num_steps + 1,
                             texp,
                             payoff_type,
                             exercise_type,
                             payoff_params)

        price = (price1 + price2) / 2.0

        return price

###############################################################################
