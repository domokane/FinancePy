# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""


from ..products.FinOptionTypes import FinOptionTypes

import numpy as np
from numba import njit, float64, int64

bump = 1e-4
from ...finutils.FinHelperFunctions import labelToString


###############################################################################
###############################################################################

@njit(
    float64[:](
        float64,
        float64,
        float64,
        float64,
        int64,
        float64,
        int64,
        float64,
        int64),
    fastmath=True,
    cache=True)
def crrTreeVal(stockPrice,
               ccInterestRate,  # continuously compounded
               ccDividendRate,  # continuously compounded
               volatility,  # Black scholes volatility
               numStepsPerYear,
               timeToExpiry,
               optionType,
               strikePrice,
               isEven):
    ''' Value an American option using a Binomial Treee '''

    numSteps = int(numStepsPerYear * timeToExpiry)

    if numSteps < 30:
        numSteps = 30

    ## OVERRIDE JUST TO SEE
    numSteps = numStepsPerYear

    # if the number of steps is even but we want odd then make it odd
    if numSteps % 2 == 0 and isEven == 0:
        numSteps += 1
    elif numSteps % 2 == 1 and isEven == 1:
        numSteps += 1

#    print(numSteps)
    # this is the size of the step
    dt = timeToExpiry / numSteps
    r = ccInterestRate
    q = ccDividendRate

    # the number of nodes on the tree
    numNodes = int(0.5 * (numSteps + 1) * (numSteps + 2))
    stockValues = np.zeros(numNodes)
    stockValues[0] = stockPrice

    optionValues = np.zeros(numNodes)
    u = np.exp(volatility * np.sqrt(dt))
    d = 1.0 / u
    sLow = stockPrice

    probs = np.zeros(numSteps)
    periodDiscountFactors = np.zeros(numSteps)

    # store time independent information for later use in tree
    for iTime in range(0, numSteps):
        a = np.exp((r - q) * dt)
        probs[iTime] = (a - d) / (u - d)
        periodDiscountFactors[iTime] = np.exp(-r * dt)

    for iTime in range(1, numSteps + 1):
        sLow *= d
        s = sLow
        for iNode in range(0, iTime + 1):
            index = 0.5 * iTime * (iTime + 1)
            stockValues[int(index + iNode)] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * numSteps * (numSteps + 1))

    for iNode in range(0, iTime + 1):

        s = stockValues[index + iNode]

        if optionType == FinOptionTypes.EUROPEAN_CALL.value:
            optionValues[index + iNode] = max(s - strikePrice, 0)
        elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
            optionValues[index + iNode] = max(strikePrice - s, 0)
        elif optionType == FinOptionTypes.AMERICAN_CALL.value:
            optionValues[index + iNode] = max(s - strikePrice, 0)
        elif optionType == FinOptionTypes.AMERICAN_PUT.value:
            optionValues[index + iNode] = max(strikePrice - s, 0)

    # begin backward steps from expiry to value date
    for iTime in range(numSteps - 1, -1, -1):

        index = int(0.5 * iTime * (iTime + 1))

        for iNode in range(0, iTime + 1):

            s = stockValues[index + iNode]

            exerciseValue = 0.0

            if optionType == FinOptionTypes.EUROPEAN_CALL.value:
                exerciseValue = 0.0
            elif optionType == FinOptionTypes.EUROPEAN_PUT.value:
                exerciseValue = 0.0
            elif optionType == FinOptionTypes.AMERICAN_CALL.value:
                exerciseValue = max(s - strikePrice, 0.0)
            elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                exerciseValue = max(strikePrice - s, 0.0)

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
                optionValues[index + iNode] = max(exerciseValue, holdValue)
            elif optionType == FinOptionTypes.AMERICAN_PUT.value:
                optionValues[index + iNode] = max(exerciseValue, holdValue)

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
###############################################################################

def crrTreeValAvg(stockPrice,
                  ccInterestRate,  # continuously compounded
                  ccDividendRate,  # continuously compounded
                  volatility,  # Black scholes volatility
                  numStepsPerYear,
                  timeToExpiry,
                  optionType,
                  strikePrice):

        value1 = crrTreeVal(stockPrice,
                            ccInterestRate,
                            ccDividendRate,
                            volatility,
                            numStepsPerYear,
                            timeToExpiry,
                            optionType,
                            strikePrice,
                            1)  # even

        value2 = crrTreeVal(stockPrice,
                            ccInterestRate,
                            ccDividendRate,
                            volatility,
                            numStepsPerYear,
                            timeToExpiry,
                            optionType,
                            strikePrice,
                            0)  # odd

        v = (value1 + value2) / 2.0
        res = {'value': v[0], 'delta': v[1], 'gamma': v[2], 'theta': v[3]}
        return res

###############################################################################
###############################################################################
