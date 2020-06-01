# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.products.equity.FinEquityBasketOption import FinEquityBasketOption
from financepy.finutils.FinOptionTypes import FinOptionTypes
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.finutils.FinDate import FinDate
import numpy as np
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBasketOption():

    import time

    valueDate = FinDate(2015, 1, 1)
    expiryDate = FinDate(2016, 1, 1)
    volatility = 0.30
    interestRate = 0.05
    discountCurve = FinFlatCurve(valueDate, interestRate)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.ones(numAssets) * volatility
    dividendYields = np.ones(numAssets) * 0.01
    stockPrices = np.ones(numAssets) * 100

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for numPaths in [10000]:
            callOption = FinEquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, numAssets)
            betas = np.ones(numAssets) * beta
            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas)
            vMC = callOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividendYields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stockPrices = np.array([100, 105, 120, 100, 90])

    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for numPaths in [10000]:

            callOption = FinEquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, numAssets)
            betas = np.ones(numAssets) * beta
            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas)
            vMC = callOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # Homogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.ones(numAssets) * volatility
    dividendYields = np.ones(numAssets) * 0.01
    stockPrices = np.ones(numAssets) * 100
    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:
        for numPaths in [10000]:
            callOption = FinEquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT, numAssets)
            betas = np.ones(numAssets) * beta
            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas)
            vMC = callOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)

    ##########################################################################
    # INHomogeneous Basket
    ##########################################################################

    numAssets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividendYields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stockPrices = np.array([100, 105, 120, 100, 90])
    betaList = np.linspace(0.0, 0.999999, 11)

    testCases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in betaList:

        for numPaths in [10000]:

            callOption = FinEquityBasketOption(
                expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT, numAssets)
            betas = np.ones(numAssets) * beta
            start = time.time()
            v = callOption.value(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas)
            vMC = callOption.valueMC(
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                betas,
                numPaths)
            end = time.time()
            duration = end - start
            testCases.print(numPaths, beta, v, vMC, duration)


test_FinBasketOption()
testCases.compareTestCases()
