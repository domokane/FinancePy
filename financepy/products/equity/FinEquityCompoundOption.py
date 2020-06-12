##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt

import numpy as np
from scipy import optimize
from numba import njit
from typing import Union

from ...finutils.FinMath import N, phi2
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...finutils.FinOptionTypes import FinOptionTypes

from ...products.equity.FinEquityOption import FinEquityOption
from ...products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from ...market.curves.FinFlatCurve import FinFlatCurve
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes 
from ...finutils.FinDate import FinDate

###############################################################################


def f(s0, *args):

    self = args[0]
    valueDate = args[1]
    discountCurve = args[2]
    dividendYield = args[3]
    model = args[4]
    value = args[5]

    if s0 <= 0.0:
        raise FinError("Unable to solve for stock price that fits K1")

    objFn = self.value(
        valueDate,
        s0,
        discountCurve,
        dividendYield,
        model)['value'] - value

    return objFn

###############################################################################


@njit(fastmath=True, cache=True, nogil=True)
def valueOnce(stockPrice,
              riskFreeRate,
              dividendYield,
              volatility,
              t1,
              t2,
              optionType1,
              optionType2,
              k1,
              k2,
              numSteps):

    if numSteps < 3:
        numSteps = 3

    # Need equally spaced time intervals for a recombining tree
    # Downside is that we may not measure periods exactly
    dt = t2 / numSteps
    numSteps1 = int(t1 / dt)
    numSteps2 = numSteps - numSteps1
    dt1 = dt
    dt2 = dt
    r = riskFreeRate
    q = dividendYield

    # print("T1:",t1,"T2:",t2,"dt:",dt,"N1*dt",numSteps1*dt,"N*dt",numSteps*dt)

    # the number of nodes on the tree
    numNodes = int(0.5 * (numSteps + 1) * (numSteps + 2))
    optionValues = np.zeros(numNodes)

    u1 = np.exp(volatility * np.sqrt(dt))
    d1 = 1.0 / u1
    u2 = np.exp(volatility * np.sqrt(dt))
    d2 = 1.0 / u2

    probs = np.zeros(numSteps)
    periodDiscountFactors = np.zeros(numSteps)

    # store time independent information for later use in tree
    for iTime in range(0, numSteps1):
        a1 = np.exp((r - q) * dt1)
        probs[iTime] = (a1 - d1) / (u1 - d1)
        periodDiscountFactors[iTime] = np.exp(-r * dt1)

    for iTime in range(numSteps1, numSteps):
        a2 = np.exp((r - q) * dt2)
        probs[iTime] = (a2 - d2) / (u2 - d2)
        periodDiscountFactors[iTime] = np.exp(-r * dt2)

    stockValues = np.zeros(numNodes)
    stockValues[0] = stockPrice
    sLow = stockPrice

    for iTime in range(1, numSteps1 + 1):
        sLow *= d1
        s = sLow
        for iNode in range(0, iTime + 1):
            index = int(0.5 * iTime * (iTime + 1))
            stockValues[index + iNode] = s
            s = s * (u1 * u1)

    for iTime in range(numSteps1 + 1, numSteps + 1):
        sLow *= d2
        s = sLow
        for iNode in range(0, iTime + 1):
            index = int(0.5 * iTime * (iTime + 1))
            stockValues[index + iNode] = s
            s = s * (u2 * u2)

    # work backwards by first setting values at expiry date t2
    index = int(0.5 * numSteps * (numSteps + 1))

    for iNode in range(0, iTime + 1):
        s = stockValues[index + iNode]
        if optionType2 == FinOptionTypes.EUROPEAN_CALL:
            optionValues[index + iNode] = max(s - k2, 0)
        elif optionType2 == FinOptionTypes.EUROPEAN_PUT:
            optionValues[index + iNode] = max(k2 - s, 0)

    # begin backward steps from expiry at t2 to first expiry at time t1
    for iTime in range(numSteps - 1, numSteps1, -1):
        index = int(0.5 * iTime * (iTime + 1))
        for iNode in range(0, iTime + 1):
            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1
            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue
            optionValues[index + iNode] = holdValue

    # Now do payoff at the end of the first expiry period at t1
    iTime = numSteps1
    index = int(0.5 * iTime * (iTime + 1))
    for iNode in range(0, iTime + 1):
        nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
        nextNodeDn = nextIndex + iNode
        nextNodeUp = nextIndex + iNode + 1
        vUp = optionValues[nextNodeUp]
        vDn = optionValues[nextNodeDn]
        futureExpectedValue = probs[iTime] * vUp
        futureExpectedValue += (1.0 - probs[iTime]) * vDn
        holdValue = periodDiscountFactors[iTime] * futureExpectedValue

        if optionType1 == FinOptionTypes.EUROPEAN_CALL:
            optionValues[index + iNode] = max(holdValue - k1, 0)
        elif optionType1 == FinOptionTypes.EUROPEAN_PUT:
            optionValues[index + iNode] = max(k1 - holdValue, 0)

    # begin backward steps from t1 expiry to value date
    for iTime in range(numSteps1 - 1, -1, -1):
        index = int(0.5 * iTime * (iTime + 1))
        for iNode in range(0, iTime + 1):
            s = stockValues[index + iNode]
            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1
            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue
            optionValues[index + iNode] = holdValue

    verbose = False
    if verbose:
        print("numSteps1:", numSteps1)
        print("numSteps2:", numSteps2)
        print("u1:", u1, "u2:", u2)
        print("dfs", periodDiscountFactors)
        print("probs:", probs)
        print("s:", stockValues)
        print("v4:", optionValues)

    # We calculate all of the important Greeks in one go
    price = optionValues[0]
    delta = (optionValues[2] - optionValues[1]) / \
        (stockValues[2] - stockValues[1])
    deltaUp = (optionValues[5] - optionValues[4]) / \
        (stockValues[5] - stockValues[4])
    deltaDn = (optionValues[4] - optionValues[3]) / \
        (stockValues[4] - stockValues[3])
    gamma = (deltaUp - deltaDn) / (stockValues[2] - stockValues[1])
    theta = (optionValues[4] - optionValues[0]) / (2.0 * dt1)
    results = np.array([price, delta, gamma, theta])
    return results

###############################################################################


class FinEquityCompoundOption(FinEquityOption):

    def __init__(self,
                 expiryDate1: FinDate,
                 expiryDate2: FinDate,
                 strikePrice1: Union[int, float],
                 strikePrice2: Union[int, float],
                 optionType1: FinOptionTypes,
                 optionType2: FinOptionTypes):

        checkArgumentTypes(self.__init__, locals())

        if expiryDate1 > expiryDate2:
            raise FinError("Expiry date 1 must preced expiry date 2")

        self._expiryDate1 = expiryDate1
        self._expiryDate2 = expiryDate2
        self._strikePrice1 = float(strikePrice1)
        self._strikePrice2 = float(strikePrice2)
        self._optionType1 = optionType1
        self._optionType2 = optionType2

    def value(self,
              valueDate,
              stockPrice,
              discountCurve,
              dividendYield,
              model):

        t1 = (self._expiryDate1 - valueDate) / gDaysInYear
        t2 = (self._expiryDate2 - valueDate) / gDaysInYear

        df2 = discountCurve.df(t2)
        r = -log(df2)/t2

        if abs(t1) < gSmall:
            t1 = gSmall

        if abs(t2) < gSmall:
            t2 = gSmall

        v = model._volatility

        if abs(v) < gSmall:
            v = gSmall

        s0 = stockPrice
        q = dividendYield
        k1 = self._strikePrice1
        k2 = self._strikePrice2

        sstar = self.impliedStockPrice(
            s0,
            self._expiryDate1,
            self._expiryDate2,
            k1,
            k2,
            self._optionType2,
            r,
            q,
            model)

        a1 = (log(s0 / sstar) + (r - q + (v**2) / 2.0) * t1) / v / sqrt(t1)
        a2 = a1 - v * sqrt(t1)
        b1 = (log(s0 / k2) + (r - q + (v**2) / 2.0) * t2) / v / sqrt(t2)
        b2 = b1 - v * sqrt(t2)

        dq2 = exp(-q * t2)
        df1 = exp(-r * t1)
        df2 = exp(-r * t2)
        c = sqrt(t1 / t2)

        # Taken from Hull Page 532 (6th edition)

        CALL = FinOptionTypes.EUROPEAN_CALL
        PUT = FinOptionTypes.EUROPEAN_PUT

        if self._optionType1 == CALL and self._optionType2 == CALL:
            v = s0 * dq2 * phi2(a1, b1, c) - k2 * df2 * \
                phi2(a2, b2, c) - df1 * k1 * N(a2)
        elif self._optionType1 == PUT and self._optionType2 == CALL:
            v = k2 * df2 * phi2(-a2, b2, -c) - s0 * dq2 * \
                phi2(-a1, b1, -c) + df1 * k1 * N(-a2)
        elif self._optionType1 == CALL and self._optionType2 == PUT:
            v = k2 * df2 * phi2(-a2, -b2, c) - s0 * dq2 * \
                phi2(-a1, -b1, c) - df1 * k1 * N(-a2)
        elif self._optionType1 == PUT and self._optionType2 == PUT:
            v = s0 * dq2 * phi2(a1, -b1, -c) - k2 * df2 * \
                phi2(a2, -b2, -c) + df1 * k1 * N(a2)
        else:
            raise FinError("Unknown option type")

        return v

###############################################################################

    def valueTree(self,
                  valueDate,
                  stockPrice,
                  discountCurve,
                  dividendYield,
                  model,
                  numSteps=200):

        if valueDate > self._expiryDate1:
            raise FinError("Value date is after expiry date.")

        t1 = (self._expiryDate1 - valueDate) / gDaysInYear
        t2 = (self._expiryDate2 - valueDate) / gDaysInYear

        df2 = discountCurve.df(t2)
        riskFreeRate = -log(df2)/t2

        volatility = model._volatility

        v1 = valueOnce(stockPrice,
                       riskFreeRate,
                       dividendYield,
                       volatility,
                       t1, t2,
                       self._optionType1,
                       self._optionType2,
                       self._strikePrice1,
                       self._strikePrice2,
                       numSteps)

        v2 = valueOnce(stockPrice,
                       riskFreeRate,
                       dividendYield,
                       volatility,
                       t1, t2,
                       self._optionType1,
                       self._optionType2,
                       self._strikePrice1,
                       self._strikePrice2,
                       numSteps + 1)

        return (v1 + v2) / 2.0

###############################################################################

    def impliedStockPrice(self,
                          stockPrice,
                          expiryDate1,
                          expiryDate2,
                          strikePrice1,
                          strikePrice2,
                          optionType2,
                          interestRate,
                          dividendYield,
                          model):

        option = FinEquityVanillaOption(expiryDate2, strikePrice2, optionType2)

        discountCurve = FinFlatCurve(expiryDate1, interestRate)

        argtuple = (
            option,
            expiryDate1,
            discountCurve,
            dividendYield,
            model,
            strikePrice1)

        sigma = optimize.newton(
            f,
            x0=stockPrice,
            args=argtuple,
            tol=1e-8,
            maxiter=50,
            fprime2=None)
        return sigma

###############################################################################
