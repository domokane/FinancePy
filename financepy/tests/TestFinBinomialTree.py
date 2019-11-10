# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.equities.FinBinomialTree import FinBinomialTree
from financepy.products.equities.FinBinomialTree import FinTreeExerciseTypes
from financepy.products.equities.FinBinomialTree import FinTreePayoffTypes
from financepy.products.equities.FinVanillaOption import FinVanillaOption
from financepy.products.equities.FinOption import FinOptionTypes
from financepy.finutils.FinDate import FinDate

from financepy.products.equities.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve

import numpy as np
import time
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinBinomialTree():

    stockPrice = 50.0
    riskFreeRate = 0.06
    dividendYield = 0.04
    volatility = 0.40

    valueDate = FinDate(2016, 1, 1)
    expiryDate = FinDate(2017, 1, 1)

    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinFlatCurve(valueDate, riskFreeRate)

    numStepsList = [100, 500, 1000, 2000, 5000]

    strikePrice = 50.0

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = FinVanillaOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    delta = putOption.delta(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    gamma = putOption.gamma(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    theta = putOption.theta(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")

    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
        results = tree.value(
            stockPrice,
            discountCurve,
            dividendYield,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)
        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.AMERICAN
    params = np.array([-1, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")

    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()
        results = tree.value(
            stockPrice,
            discountCurve,
            dividendYield,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)
        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = FinVanillaOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    delta = callOption.delta(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    gamma = callOption.gamma(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    theta = callOption.theta(
        valueDate,
        stockPrice,
        discountCurve,
        dividendYield,
        model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()

        results = tree.value(
            stockPrice,
            discountCurve,
            dividendYield,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)

        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")

    payoff = FinTreePayoffTypes.VANILLA_OPTION
    exercise = FinTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = FinBinomialTree()

        results = tree.value(
            stockPrice,
            discountCurve,
            dividendYield,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)

        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)


test_FinBinomialTree()
testCases.compareTestCases()
