# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""

import numpy as np
from numba import njit, float64, int64

from enum import Enum

from ...finutils.FinError import FinError
from ...finutils.FinGlobalVariables import gDaysInYear
from ...products.FinOptionTypes import FinOptionTypes
from ...products.equity.FinEquityModelTypes import FinEquityModelBlackScholes

bump = 1e-4

###############################################################################


class FinImplementations(Enum):
    CRR_TREE = 1,
    BARONE_ADESI_APPROX = 2

###############################################################################
# TODO
#   Add term structure of rates inside code
#   Implement some fast approximations like Barone Adesi
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
               riskFreeRate,
               dividendYield,
               volatility,
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
    r = riskFreeRate
    q = dividendYield

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


class FinEquityAmericanOption():
    ''' Class that performs the valuation of an American style option on a
    dividend paying stock. Can easily be extended to price American style FX
    options. '''

###############################################################################

    def __init__(self,
                 expiryDate,
                 strikePrice,
                 optionType):
        ''' Create the American option using its expiry date, strike price
        and option type (call or put). The option notional is provided too. '''

        if strikePrice < 0.0:
            raise FinError("Strike price must be positive.")

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType

###############################################################################

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model,
              numStepsPerYear=100):
        ''' Value the American option using a CRR Binomial Tree. The model
        assumes that the dividend is a continuous one. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date is after expiry date.")

        # do some validation
        timeToExpiry = (self._expiryDate - valueDate) / gDaysInYear

        # Convergence is improved if I take the average of an even and
        # odd number of steps in the Binomial Tree

        volatility = model._volatility
        df = discountCurve.df(self._expiryDate)
        riskFreeRate = -np.log(df)/timeToExpiry

        value1 = crrTreeVal(stockPrice,
                            riskFreeRate,
                            dividendYield,
                            volatility,
                            numStepsPerYear,
                            timeToExpiry,
                            self._optionType.value,
                            self._strikePrice,
                            1)  # even

        value2 = crrTreeVal(stockPrice,
                            riskFreeRate,
                            dividendYield,
                            volatility,
                            numStepsPerYear,
                            timeToExpiry,
                            self._optionType.value,
                            self._strikePrice,
                            0)  # odd

        v = (value1 + value2) / 2.0
        res = {'value': v[0], 'delta': v[1], 'gamma': v[2], 'theta': v[3]}
        return res

###############################################################################

    def delta(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model):
        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice + bump,
            discountCurve,
            dividendYield,
            model)

        delta = (vBumped['value'] - v['value']) / bump
        return delta

###############################################################################

    def gamma(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model):

        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumpedDn = self.value(
            valueDate,
            stockPrice - bump,
            discountCurve,
            dividendYield,
            model)
        vBumpedUp = self.value(
            valueDate,
            stockPrice + bump,
            discountCurve,
            dividendYield,
            model)

        d2 = vBumpedUp['value'] - v['value']
        d1 = v['value'] - vBumpedDn['value']
        gamma = (d2 - d1) / bump
        return gamma

###############################################################################

    def vega(self,
             valueDate,
             stockPrice,
             discountCurve,
             dividendYield,
             model):
        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            FinEquityModelBlackScholes(model._volatility + bump))

        if type(v) is dict:
            vega = (vBumped['value'] - v['value']) / bump
        else:
            vega = (vBumped - v) / bump

        return vega

###############################################################################

    def theta(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model):

        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        nextDate = valueDate.addDays(1)
        bump = 1.0 / gDaysInYear
        vBumped = self.value(
            nextDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)

        if type(v) is dict:
            theta = (vBumped['value'] - v['value']) / bump
        else:
            theta = (vBumped - v) / bump

        return theta

###############################################################################

    def rho(self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):

        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice,
            discountCurve.bump(bump),
            dividendYield,
            model)

        if type(v) is dict:
            rho = (vBumped['value'] - v['value']) / bump
        else:
            rho = (vBumped - v) / bump

        return rho

###############################################################################
