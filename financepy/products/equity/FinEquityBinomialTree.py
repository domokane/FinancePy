##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
from enum import Enum

import numpy as np
from numba import jit, njit, float64, int64

from ...utils.FinError import FinError
from ...utils.global_variables import gDaysInYear
from ...utils.fin_math import heaviside
from ...utils.helper_functions import labelToString

###############################################################################


class FinEquityTreePayoffTypes(Enum):
    FWD_CONTRACT = 1
    VANILLA_OPTION = 2
    DIGITAL_OPTION = 3
    POWER_CONTRACT = 4
    POWER_OPTION = 5
    LOG_CONTRACT = 6
    LOG_OPTION = 7


class FinEquityTreeExerciseTypes(Enum):
    EUROPEAN = 1
    AMERICAN = 2

###############################################################################


@jit
def _validatePayoff(payoffType, payoffParams):

    numParams = 0

    if payoffType == FinEquityTreePayoffTypes.FWD_CONTRACT.value:
        numParams = 1
    elif payoffType == FinEquityTreePayoffTypes.VANILLA_OPTION.value:
        numParams = 2
    elif payoffType == FinEquityTreePayoffTypes.DIGITAL_OPTION.value:
        numParams = 2
    elif payoffType == FinEquityTreePayoffTypes.POWER_CONTRACT.value:
        numParams = 2
    elif payoffType == FinEquityTreePayoffTypes.POWER_OPTION.value:
        numParams = 3
    elif payoffType == FinEquityTreePayoffTypes.LOG_CONTRACT.value:
        numParams = 0
    elif payoffType == FinEquityTreePayoffTypes.LOG_OPTION.value:
        numParams = 1
    else:
        raise FinError("Unknown payoff type")

    if len(payoffParams) != numParams:
        raise FinError(
            "Number of parameters required for " +
            str(payoffType) +
            " must be " +
            str(numParams))

    return None

###############################################################################


@njit(float64(float64, int64, float64[:]), fastmath=True, cache=True)
def _payoffValue(s, payoffType, payoffParams):

    if payoffType == FinEquityTreePayoffTypes.FWD_CONTRACT.value:
        payoff = payoffParams[0] * s
    elif payoffType == FinEquityTreePayoffTypes.VANILLA_OPTION.value:
        payoff = max(payoffParams[0] * (s - payoffParams[1]), 0.0)
    elif payoffType == FinEquityTreePayoffTypes.DIGITAL_OPTION.value:
        payoff = heaviside(payoffParams[0] * (s - payoffParams[1]))
    elif payoffType == FinEquityTreePayoffTypes.POWER_CONTRACT.value:
        payoff = payoffParams[0] * (s**payoffParams[1])
    elif payoffType == FinEquityTreePayoffTypes.POWER_OPTION.value:
        payoff = max(payoffParams[0] *
                     ((s**payoffParams[2]) -
                      payoffParams[1]), 0.0)
    elif payoffType == FinEquityTreePayoffTypes.LOG_CONTRACT.value:
        payoff = log(s)
    elif payoffType == FinEquityTreePayoffTypes.LOG_OPTION.value:
        payoff = max(log(s) - payoffParams[0], 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff

###############################################################################


@njit(fastmath=True, cache=True)
def _valueOnce(stock_price,
               r,
               q,
               volatility,
               num_steps,
               timeToExpiry,
               payoffType,
               exerciseType,
               payoffParams):

    if num_steps < 3:
        num_steps = 3

#        validatePayoff(payoffType.value,payoffParams)

    payoffTypeValue = payoffType.value

    # this is the size of the step
    dt = timeToExpiry / num_steps

    # the number of nodes on the tree
    numNodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stockValues = np.zeros(numNodes)
    stockValues[0] = stock_price

    optionValues = np.zeros(numNodes)
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
            stockValues[int(index + iNode)] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * num_steps * (num_steps + 1))

    for iNode in range(0, iTime + 1):
        s = stockValues[index + iNode]
        optionValues[index +
                     iNode] = _payoffValue(s, payoffTypeValue, payoffParams)

    # begin backward steps from expiry
    for iTime in range(num_steps - 1, -1, -1):

        index = int(0.5 * iTime * (iTime + 1))

        for iNode in range(0, iTime + 1):

            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextNodeDn + 1
            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            if exerciseType == FinEquityTreeExerciseTypes.EUROPEAN:
                optionValues[index + iNode] = holdValue
            elif exerciseType == FinEquityTreeExerciseTypes.AMERICAN:
                s = stockValues[index + iNode]
                exerciseValue = _payoffValue(s, payoffTypeValue, payoffParams)
                optionValues[index + iNode] = max(exerciseValue, holdValue)

    price = optionValues[0]
    delta = (optionValues[2] - optionValues[1]) / \
        (stockValues[2] - stockValues[1])
    deltaUp = (optionValues[5] - optionValues[4]) / \
        (stockValues[5] - stockValues[4])
    deltaDn = (optionValues[4] - optionValues[3]) / \
        (stockValues[4] - stockValues[3])
    gamma = (deltaUp - deltaDn) / (stockValues[2] - stockValues[1])
    theta = (optionValues[4] - optionValues[0]) / (2.0 * dt)
    results = np.array([price, delta, gamma, theta])
    return results

###############################################################################


class FinEquityBinomialTree():

    def __init__(self):
        pass
#        self.m_optionValues = np.zeros()
#        self.m_stockValues = np.zeros()
#        self.m_upProbabilities = np.zeros()
#
#       self.m_num_steps = 10
#        self.m_numNodes = 10

###############################################################################

    def value(self,
              stock_price,
              discount_curve,
              dividendCurve,
              volatility,
              num_steps,
              valuation_date,
              payoff,
              expiry_date,
              payoffType,
              exerciseType,
              payoffParams):

        # do some validation
        texp = (expiry_date - valuation_date) / gDaysInYear
        r = discount_curve.zeroRate(expiry_date)

        dq = dividendCurve.df(expiry_date)
        q = -np.log(dq)/texp

        price1 = _valueOnce(stock_price,
                            r,
                            q,
                            volatility,
                            num_steps,
                            texp,
                            payoffType,
                            exerciseType,
                            payoffParams)

        # Can I reuse the same tree ?
        price2 = _valueOnce(stock_price,
                            r,
                            q,
                            volatility,
                            num_steps + 1,
                            texp,
                            payoffType,
                            exerciseType,
                            payoffParams)

        price = (price1 + price2) / 2.0

        return price

###############################################################################
