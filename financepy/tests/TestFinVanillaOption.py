# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.equities.FinOption import FinOptionTypes
from financepy.products.equities.FinVanillaOption import FinVanillaOption, FinOptionModelTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinVanillaOption():

    valueDate = FinDate(2015, 1, 1)
    expiryDate = FinDate(2015, 7, 1)
    stockPrice = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    modelType = FinOptionModelTypes.BLACKSCHOLES

    import time

    numPathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for numPaths in numPathsList:

        callOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        modelParams = (volatility)
        value = callOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        start = time.time()
        valueMC = callOption.valueMC(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)

##########################################################################

    stockPrices = range(80, 120, 2)
    numPaths = 100000

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for stockPrice in stockPrices:

        callOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        start = time.time()
        valueMC = callOption.valueMC(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)

##########################################################################

    stockPrices = range(80, 120, 2)
    numPaths = 100000

    testCases.header("STOCK PRICE", "VALUE_BS", "VALUE_MC", "TIME")

    for stockPrice in stockPrices:

        putOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT)
        value = putOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        start = time.time()
        valueMC = putOption.valueMC(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(stockPrice, value, valueMC, duration)

##########################################################################

    stockPrices = range(80, 120, 2)

    testCases.header(
        "STOCK PRICE",
        "VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for stockPrice in stockPrices:
        callOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        delta = callOption.delta(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        vega = callOption.vega(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        theta = callOption.theta(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        # callOption.rho(valueDate,stockPrice, interestRate, dividendYield, modelType, modelParams)
        rho = 999
        testCases.print(stockPrice, value, delta, vega, theta, rho)

    testCases.header(
        "STOCK PRICE",
        "VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for stockPrice in stockPrices:
        putOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT)
        value = putOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        delta = putOption.delta(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        vega = putOption.vega(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        theta = putOption.theta(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        # putOption.rho(valueDate,stockPrice, interestRate, dividendYield, modelType, modelParams)
        rho = 999
        testCases.print(stockPrice, value, delta, vega, theta, rho)

##########################################################################

    testCases.header("STOCK PRICE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    stockPrices = range(60, 150, 2)

    for stockPrice in stockPrices:
        callOption = FinVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            modelType,
            modelParams)
        impliedVol = callOption.impliedVolatility(
            valueDate, stockPrice, interestRate, dividendYield, value)
        testCases.print(stockPrice, value, volatility, impliedVol)


test_FinVanillaOption()
testCases.compareTestCases()
