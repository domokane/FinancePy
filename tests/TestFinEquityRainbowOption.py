###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from math import sqrt
import time

import sys
sys.path.append("..")

from financepy.products.equity.equity_rainbow_option import EquityRainbowOption
from financepy.products.equity.equity_rainbow_option import EquityRainbowOptionTypes
from financepy.utils.helpers import betaVectorToCorrMatrix
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityRainbowOption():

    #        import matplotlib.pyplot as plt

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    interestRate = 0.05

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    num_assets = 2
    volatilities = np.ones(num_assets) * 0.3

    dividendYields = np.ones(num_assets) * 0.01

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

    stockPrices = np.ones(num_assets) * 100
    numPathsList = [10000]
    corrList = np.linspace(0.0, 0.999999, 6)
    strike = 100.0

    testCases.banner(
        "===================================================================")
    testCases.banner("                      CALL ON MAXIMUM")
    testCases.banner(
        "===================================================================")

    payoff_type = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiryDate, payoff_type, payoff_params, num_assets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()
            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiryDate, payoff_type, payoff_params, num_assets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()

            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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

    payoff_type = EquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiryDate, payoff_type, payoff_params, num_assets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()

            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiryDate, payoff_type, payoff_params, num_assets)

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    testCases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            start = time.time()
            v = rainbowOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)
            v_MC = rainbowOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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

    num_assets = 2
    volatilities = np.ones(num_assets) * 0.3
    dividendYields = np.ones(num_assets) * 0.01
    stockPrices = np.ones(num_assets) * 100
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

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
            payoff_params1 = [strike]
            rainbowOption1 = EquityRainbowOption(
                expiryDate, payoff_type1, payoff_params1, num_assets)

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [1, strike]
            rainbowOption2 = EquityRainbowOption(
                expiryDate, payoff_type2, payoff_params2, num_assets)

            start = time.time()

            v = rainbowOption1.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption2.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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

        betas = np.ones(num_assets) * sqrt(correlation)
        corrMatrix = betaVectorToCorrMatrix(betas)

        for numPaths in numPathsList:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MINIMUM
            payoff_params1 = [strike]
            rainbowOption1 = EquityRainbowOption(
                expiryDate, payoff_type1, payoff_params1, num_assets)

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [2, strike]
            rainbowOption2 = EquityRainbowOption(
                expiryDate, payoff_type2, payoff_params2, num_assets)

            start = time.time()

            v = rainbowOption1.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            v_MC = rainbowOption2.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividendYields = np.ones(num_assets) * 0.01
    stockPrices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

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

        payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
        payoff_params2 = [n, strike]
        rainbowOption2 = EquityRainbowOption(
            expiryDate, payoff_type2, payoff_params2, num_assets)

        for correlation in corrList:

            betas = np.ones(num_assets) * sqrt(correlation)
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v_MC = rainbowOption2.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividendYields = np.ones(num_assets) * 0.01
    stockPrices = np.ones(num_assets) * 100

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

        payoff_type2 = EquityRainbowOptionTypes.PUT_ON_NTH
        payoff_params2 = [n, strike]
        rainbowOption2 = EquityRainbowOption(
            expiryDate, payoff_type2, payoff_params2, num_assets)

        for correlation in corrList:

            betas = np.ones(num_assets) * sqrt(correlation)
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v_MC = rainbowOption2.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
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


test_EquityRainbowOption()
testCases.compareTestCases()
