###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.products.equity.equity_basket_option import EquityBasketOption
from financepy.utils.global_types import FinOptionTypes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.helpers import betaVectorToCorrMatrix
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityBasketOption():

    import time

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    volatility = 0.30
    interestRate = 0.05
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividendYields = np.ones(num_assets) * 0.01
    stockPrices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for numPaths in [10000]:
            callOption = EquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)

            vMC = callOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividendYields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stockPrices = np.array([100, 105, 120, 100, 90])

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for numPaths in [10000]:

            callOption = EquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v = callOption.value(
                    valueDate,
                    stockPrices,
                    discountCurve,
                    dividend_curves,
                    volatilities,
                    corrMatrix)

            vMC = callOption.value_mc(
                    valueDate,
                    stockPrices,
                    discountCurve,
                    dividend_curves,
                    volatilities,
                    corrMatrix,
                    numPaths)

            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividendYields = np.ones(num_assets) * 0.01
    stockPrices = np.ones(num_assets) * 100
    betaList = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for numPaths in [10000]:
            callOption = EquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)
            vMC = callOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividendYields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stockPrices = np.array([100, 105, 120, 100, 90])
    betaList = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividendYields:
        dividend_curve = FinDiscountCurveFlat(valueDate, q)
        dividend_curves.append(dividend_curve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for numPaths in [10000]:

            callOption = EquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix)
            vMC = callOption.value_mc(
                valueDate,
                stockPrices,
                discountCurve,
                dividend_curves,
                volatilities,
                corrMatrix,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)


###############################################################################


test_EquityBasketOption()
testCases.compareTestCases()
