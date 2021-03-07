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

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    volatility = 0.30
    interestRate = 0.05
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for num_paths in [10000]:
            callOption = EquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix)

            vMC = callOption.value_mc(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for num_paths in [10000]:

            callOption = EquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()

            v = callOption.value(
                    valuation_date,
                    stock_prices,
                    discount_curve,
                    dividend_curves,
                    volatilities,
                    corrMatrix)

            vMC = callOption.value_mc(
                    valuation_date,
                    stock_prices,
                    discount_curve,
                    dividend_curves,
                    volatilities,
                    corrMatrix,
                    num_paths)

            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100
    betaList = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for num_paths in [10000]:
            callOption = EquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix)
            vMC = callOption.value_mc(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])
    betaList = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for num_paths in [10000]:

            callOption = EquityBasketOption(
                expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, num_assets)
            betas = np.ones(num_assets) * beta
            corrMatrix = betaVectorToCorrMatrix(betas)

            start = time.time()
            v = callOption.value(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix)
            vMC = callOption.value_mc(
                valuation_date,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corrMatrix,
                num_paths)
            end = time.time()
            duration = end - start
            testCases.print(num_paths, beta, v, vMC, duration)


###############################################################################


test_EquityBasketOption()
testCases.compareTestCases()
