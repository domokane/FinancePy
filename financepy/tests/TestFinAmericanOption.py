# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:04:57 2019

@author: Dominic
"""

from financepy.products.equities.FinVanillaOption import FinVanillaOption
from financepy.products.equities.FinVanillaOption import FinOptionTypes
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.equities.FinAmericanOption import FinAmericanOption
from financepy.products.equities.FinVanillaOption import FinOptionModelTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def testFinAmericanOption():

    import time

    stockPrice = 50.0
    riskFreeRate = 0.06
    dividendYield = 0.04
    volatility = 0.40

    strikePrice = 50.0
    valueDate = FinDate(2016, 1, 1)
    expiryDate = FinDate(2017, 1, 1)

    numStepsList = [1000, 2000]

    testCases.banner("================== EUROPEAN PUT =======================")
    modelType = FinOptionModelTypes.BLACKSCHOLES
    modelParams = (volatility)
    putOption = FinVanillaOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    delta = putOption.delta(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    gamma = putOption.gamma(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    theta = putOption.theta(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = FinAmericanOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for numSteps in numStepsList:
        start = time.time()
        results = option.value(
            valueDate,
            stockPrice,
            riskFreeRate,
            dividendYield,
            volatility,
            numSteps)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", numSteps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = FinAmericanOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for numSteps in numStepsList:
        start = time.time()
        results = option.value(
            valueDate,
            stockPrice,
            riskFreeRate,
            dividendYield,
            volatility,
            numSteps)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", numSteps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = FinVanillaOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    delta = callOption.delta(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    gamma = callOption.gamma(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)
    theta = callOption.theta(
        valueDate,
        stockPrice,
        riskFreeRate,
        dividendYield,
        modelType,
        modelParams)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = FinAmericanOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for numSteps in numStepsList:
        start = time.time()
        results = option.value(
            valueDate,
            stockPrice,
            riskFreeRate,
            dividendYield,
            volatility,
            numSteps)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_CALL", numSteps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")
    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    option = FinAmericanOption(
        expiryDate,
        strikePrice,
        FinOptionTypes.AMERICAN_CALL)

    for numSteps in numStepsList:
        start = time.time()
        results = option.value(
            valueDate,
            stockPrice,
            riskFreeRate,
            dividendYield,
            volatility,
            numSteps)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_CALL", numSteps, results, duration)

#    FinTest.TestReport(filename)


testFinAmericanOption()
testCases.compareTestCases()
