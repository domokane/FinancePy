##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ..utils.global_types import FinOptionTypes

import numpy as np
from numba import njit, float64, int64

bump = 1e-4

###############################################################################
###############################################################################

@njit(float64[:](float64, float64, float64, float64, int64, float64, int64,
                 float64, int64), fastmath=True, cache=True)
def crrTreeVal(stock_price,
               ccInterestRate,  # continuously compounded
               ccDividendRate,  # continuously compounded
               volatility,  # Black scholes volatility
               num_steps_per_year,
               timeToExpiry,
               optionType,
               strikePrice,
               isEven):
    """ Value an American option using a Binomial Treee """

    num_steps = int(num_steps_per_year * timeToExpiry)

    if num_steps < 30:
        num_steps = 30

    ## OVERRIDE JUST TO SEE
    num_steps = num_steps_per_year

    # if the number of steps is even but we want odd then make it odd
    if num_steps % 2 == 0 and isEven == 0:
        num_steps += 1
    elif num_steps % 2 == 1 and isEven == 1:
        num_steps += 1

#    print(num_steps)
    # this is the size of the step
    dt = timeToExpiry / num_steps
    r = ccInterestRate
    q = ccDividendRate

    # the number of nodes on the tree
    numNodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stockValues = np.zeros(numNodes)
    stockValues[0] = stock_price

    optionValues = np.zeros(numNodes)
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
            stockValues[int(index + iNode)] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * num_steps * (num_steps + 1))

    for iNode in range(0, iTime + 1):

        s = stockValues[index + iNode]

        if optionType == FinOptionTypes.EUROPEAN_CALL.value:
            optionValues[index + iNode] = np.maximum(s - strikePrice, 0.0)
        elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
            optionValues[index + iNode] = np.maximum(strikePrice - s, 0.0)
        elif optionType == FinOptionTypes.AMERICAN_CALL.value:
            optionValues[index + iNode] = np.maximum(s - strikePrice, 0.0)
        elif optionType == FinOptionTypes.AMERICAN_PUT.value:
            optionValues[index + iNode] = np.maximum(strikePrice - s, 0.0)

    # begin backward steps from expiry to value date
    for iTime in range(num_steps - 1, -1, -1):

        index = int(0.5 * iTime * (iTime + 1))

        for iNode in range(0, iTime + 1):

            s = stockValues[index + iNode]

            exerciseValue = 0.0

            if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                exerciseValue = 0.0
            elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                exerciseValue = 0.0
            elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                exerciseValue = np.maximum(s - strikePrice, 0.0)
            elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                exerciseValue = np.maximum(strikePrice - s, 0.0)

            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))

            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1

            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                optionValues[index + iNode] = holdValue
            elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                optionValues[index + iNode] = holdValue
            elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                optionValues[index + iNode] = np.maximum(exerciseValue, holdValue)
            elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                optionValues[index + iNode] = np.maximum(exerciseValue, holdValue)

    # We calculate all of the important Greeks in one go
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


def crrTreeValAvg(stock_price,
                  ccInterestRate,  # continuously compounded
                  ccDividendRate,  # continuously compounded
                  volatility,  # Black scholes volatility
                  num_steps_per_year,
                  timeToExpiry,
                  optionType,
                  strikePrice):
    """ Calculate the average values off the tree using an even and an odd
    number of time steps. """

    value1 = crrTreeVal(stock_price,
                        ccInterestRate,
                        ccDividendRate,
                        volatility,
                        num_steps_per_year,
                        timeToExpiry,
                        optionType,
                        strikePrice,
                        1)  # even

    value2 = crrTreeVal(stock_price,
                        ccInterestRate,
                        ccDividendRate,
                        volatility,
                        num_steps_per_year,
                        timeToExpiry,
                        optionType,
                        strikePrice,
                        0)  # odd

    v = (value1 + value2) / 2.0
    res = {'value': v[0], 'delta': v[1], 'gamma': v[2], 'theta': v[3]}
    return res

###############################################################################
