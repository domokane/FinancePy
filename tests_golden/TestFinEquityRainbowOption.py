###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from math import sqrt
import time

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityRainbowOption import FinEquityRainbowOption
from financepy.products.equity.FinEquityRainbowOption import FinEquityRainbowOptionTypes
from financepy.finutils.FinHelperFunctions import betaVectorToCorrMatrix
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityRainbowOption():

    #        import matplotlib.pyplot as plt

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    interestRate = 0.05

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    numAssets = 2
    volatilities = np.ones(numAssets) * 0.3

    dividendYields = np.ones(numAssets) * 0.01

    dividendCurves = []
    for q in dividendYields:
        dividendCurve = FinDiscountCurveFlat(valueDate, q)
        dividendCurves.append(dividendCurve)

    stockPrices = np.ones(numAssets) * 100
    numPathsList = [10000]
    corrList = np.linspace(0.0, 0.999999, 6)
    strike = 100.0

    testCases.banner(
        "===================================================================")
    testCases.banner("                      CALL ON MAXIMUM")
    testCases.banner(
        "===================================================================")

    payoffType = FinEquityRainbowOptionTypes.CALL_ON_MAXIMUM
    payoffParams = [strike]
    rainbowOption = FinEquityRainbowOption(
        expiryDate, payoffType, payoffParams, numAssets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()
            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON MAX Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner("                       CALL ON MINIMUM")
    testCases.banner(
        "===================================================================")
    payoffType = FinEquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoffParams = [strike]
    rainbowOption = FinEquityRainbowOption(
        expiryDate, payoffType, payoffParams, numAssets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()

            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON MIN Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

###############################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner("                      PUT ON MAXIMUM")
    testCases.banner(
        "===================================================================")

    payoffType = FinEquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoffParams = [strike]
    rainbowOption = FinEquityRainbowOption(
        expiryDate, payoffType, payoffParams, numAssets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()

            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "PUT ON MAX Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "PUT ON MAX Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

##########################################################################

    testCases.banner(
        "===================================================================")
    testCases.banner("                       PUT ON MINIMUM")
    testCases.banner(
        "===================================================================")
    payoffType = FinEquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoffParams = [strike]
    rainbowOption = FinEquityRainbowOption(
        expiryDate, payoffType, payoffParams, numAssets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()
            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)
            v_MC = rainbowOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "PUT ON MIN Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "PUT ON MIN Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

##########################################################################

    numAssets = 2
    volatilities = np.ones(numAssets) * 0.3
    dividendYields = np.ones(numAssets) * 0.01
    stockPrices = np.ones(numAssets) * 100
    strike = 100.0
    correlation = 0.50

    testCases.banner(
        "===================================================================")
    testCases.banner("                      CALL ON 1st")
    testCases.banner(
        "===================================================================")

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            payoffType1 = FinEquityRainbowOptionTypes.CALL_ON_MAXIMUM
            payoffParams1 = [strike]
            rainbowOption1 = FinEquityRainbowOption(
                expiryDate, payoffType1, payoffParams1, numAssets)

            payoffType2 = FinEquityRainbowOptionTypes.CALL_ON_NTH
            payoffParams2 = [1, strike]
            rainbowOption2 = FinEquityRainbowOption(
                expiryDate, payoffType2, payoffParams2, numAssets)

            start = time.time()

            v = rainbowOption1.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption2.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON 1st Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

    testCases.banner(
        "===================================================================")
    testCases.banner("                      CALL ON 2nd")
    testCases.banner(
        "===================================================================")

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(numAssets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            payoffType1 = FinEquityRainbowOptionTypes.CALL_ON_MINIMUM
            payoffParams1 = [strike]
            rainbowOption1 = FinEquityRainbowOption(
                expiryDate, payoffType1, payoffParams1, numAssets)

            payoffType2 = FinEquityRainbowOptionTypes.CALL_ON_NTH
            payoffParams2 = [2, strike]
            rainbowOption2 = FinEquityRainbowOption(
                expiryDate, payoffType2, payoffParams2, numAssets)

            start = time.time()

            v = rainbowOption1.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption2.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

#    plt.figure(figsize=(10,8))
#    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
#    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON 2nd Rainbow Option MC")
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

    testCases.banner(
        "===================================================================")
    testCases.banner("                      CALL ON 1-5")
    testCases.banner(
        "===================================================================")

    rainboxOptionValues = []
    rainbowOptionValuesMC = []
    numPaths = 10000
    numAssets = 5
    volatilities = np.ones(numAssets) * 0.3
    dividendYields = np.ones(numAssets) * 0.01
    stockPrices = np.ones(numAssets) * 100

    dividendCurves = []
    for q in dividendYields:
        dividendCurve = FinDiscountCurveFlat(valueDate, q)
        dividendCurves.append(dividendCurve)

#    plt.figure(figsize=(10,8))

    testCases.header(
        "NUMPATHS",
        "CORRELATION",
        "NTD",
        "VALUE",
        "VALUE_MC",
        "TIME")

    for n in [1, 2, 3, 4, 5]:

        rainboxOptionValues = []
        rainbowOptionValuesMC = []

        payoffType2 = FinEquityRainbowOptionTypes.CALL_ON_NTH
        payoffParams2 = [n, strike]
        rainbowOption2 = FinEquityRainbowOption(
            expiryDate, payoffType2, payoffParams2, numAssets)

        for correlation in corrList:

            betas = np.ones(numAssets) * sqrt(correlation)
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v_MC = rainbowOption2.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, n, v, v_MC, duration)

            rainbowOptionValuesMC.append(v_MC)

#        plt.plot(corrList, rainbowOptionValuesMC, 'o-', label = "CALL Rainbow Option MC NTH = " + str(n))
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')

    testCases.banner(
        "===================================================================")
    testCases.banner("                      PUT ON 1-5")
    testCases.banner(
        "===================================================================")

    rainboxOptionValues = []
    rainbowOptionValuesMC = []
    numPaths = 10000
    numAssets = 5
    volatilities = np.ones(numAssets) * 0.3
    dividendYields = np.ones(numAssets) * 0.01
    stockPrices = np.ones(numAssets) * 100

#    plt.figure(figsize=(10,8))

    testCases.header(
        "NUMPATHS",
        "CORRELATION",
        "NTD",
        "VALUE",
        "VALUE_MC",
        "TIME")

    for n in [1, 2, 3, 4, 5]:

        rainboxOptionValues = []
        rainbowOptionValuesMC = []

        payoffType2 = FinEquityRainbowOptionTypes.PUT_ON_NTH
        payoffParams2 = [n, strike]
        rainbowOption2 = FinEquityRainbowOption(
            expiryDate, payoffType2, payoffParams2, numAssets)

        for correlation in corrList:

            betas = np.ones(numAssets) * sqrt(correlation)
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v_MC = rainbowOption2.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendCurves,
                volatilities,
                corrMatrix,
                numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, correlation, n, v, v_MC, duration)

            rainbowOptionValuesMC.append(v_MC)

#    plt.plot(corrList, rainbowOptionValuesMC, 'o-', label = "PUT Rainbow Option MC NTH = " + str(n))
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')


###############################################################################


test_FinEquityRainbowOption()
testCases.compareTestCases()
