###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityBasketOption import FinEquityBasketOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.helper_functions import betaVectorToCorrMatrix
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityBasketOption():

    import time

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    volatility = 0.30
    interestRate = 0.05
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.ones(numAssets) * volatility
    dividend_yields = np.ones(numAssets) * 0.01
    stock_prices = np.ones(numAssets) * 100

    dividendCurves = []
    for q in dividend_yields:
        dividendCurve = DiscountCurveFlat(valuation_date, q)
        dividendCurves.append(dividendCurve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for num_paths in [10000]:
            callOption = FinEquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, numAssets)
            betas = np.ones(numAssets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix)

            vMC = callOption.valueMC(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])

    dividendCurves = []
    for q in dividend_yields:
        dividendCurve = DiscountCurveFlat(valuation_date, q)
        dividendCurves.append(dividendCurve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for num_paths in [10000]:

            callOption = FinEquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, numAssets)
            betas = np.ones(numAssets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v = callOption.value(
                    valuation_date,
                    stock_prices,
                    discount_curve,
                    dividendCurves,
                    volatilities,
                    corrMatrix)

            vMC = callOption.valueMC(
                    valuation_date,
                    stock_prices,
                    discount_curve,
                    dividendCurves,
                    volatilities,
                    corrMatrix,
                    num_paths)

            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.ones(numAssets) * volatility
    dividend_yields = np.ones(numAssets) * 0.01
    stock_prices = np.ones(numAssets) * 100
    betaList = np.linspace(0.0, 0.999999, 11)

    dividendCurves = []
    for q in dividend_yields:
        dividendCurve = DiscountCurveFlat(valuation_date, q)
        dividendCurves.append(dividendCurve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for num_paths in [10000]:
            callOption = FinEquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, numAssets)
            betas = np.ones(numAssets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix)
            vMC = callOption.valueMC(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])
    betaList = np.linspace(0.0, 0.999999, 11)

    dividendCurves = []
    for q in dividend_yields:
        dividendCurve = DiscountCurveFlat(valuation_date, q)
        dividendCurves.append(dividendCurve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for num_paths in [10000]:

            callOption = FinEquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, numAssets)
            betas = np.ones(numAssets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix)
            vMC = callOption.valueMC(
                valuation_date,
                stock_prices,
                discount_curve,
                dividendCurves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)


###############################################################################


test_FinEquityBasketOption()
testCases.compareTestCases()
