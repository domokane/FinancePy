# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import time
import numpy as np

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.fx.FinFXOption import FinFXOptionTypes
from financepy.products.fx.FinFXVanillaOption import FinFXVanillaOption
from financepy.products.fx.FinFXModelTypes import FinFXModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve

from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXVanillaOption():
#   Example from Hull 4th edition page 284
    valueDate = FinDate(2015, 1, 1)
    expiryDate = valueDate.addMonths(4)
    spotFXRate = 1.60
    volatility = 0.1411
    domInterestRate = 0.08
    forInterestRate = 0.11
    model = FinFXModelBlackScholes(volatility)
    domDiscountCurve = FinFlatCurve(valueDate, domInterestRate)
    forDiscountCurve = FinFlatCurve(valueDate, forInterestRate)

    numPathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")
    strikeFXRate = 1.60

    for numPaths in numPathsList:

        callOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_CALL)

        value = callOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)

        start = time.time()

        valueMC = callOption.valueMC(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model,
            numPaths)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)

##########################################################################

    spotFXRates = np.arange(100, 200, 10)
    spotFXRates = spotFXRates/100.0
    numPaths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC", "TIME")

    for spotFXRate in spotFXRates:

        callOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        start = time.time()
        valueMC = callOption.valueMC(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)

##########################################################################

    spotFXRates = np.arange(100, 200, 10) / 100.0
    numPaths = 100000

    testCases.header("SPOT FX RATE", "PUT_VALUE_BS", "PUT_VALUE_MC", "TIME")

    for spotFXRate in spotFXRates:

        putOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_PUT)
        value = putOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        start = time.time()
        valueMC = putOption.valueMC(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(spotFXRate, value, valueMC, duration)

##########################################################################

    spotFXRates = np.arange(100, 200, 10)/100.0

    testCases.header(
        "SPOT FX RATE",
        "CALL_VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for spotFXRate in spotFXRates:
        callOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        delta = callOption.delta(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        vega = callOption.vega(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        theta = callOption.theta(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        # callOption.rho(valueDate,stockPrice, interestRate, dividendYield, modelType, modelParams)
        rho = 999
        testCases.print(spotFXRate, value, delta, vega, theta, rho)

    testCases.header(
        "SPOT FX RATE",
        "PUT_VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for spotFXRate in spotFXRates:
        putOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_PUT)
        value = putOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        delta = putOption.delta(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        vega = putOption.vega(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        theta = putOption.theta(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        # putOption.rho(valueDate,stockPrice, interestRate, dividendYield, 
        # modelType, modelParams)
        rho = 999
        testCases.print(spotFXRate, value, delta, vega, theta, rho)

##########################################################################

    testCases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spotFXRates = np.arange(100, 200, 10)/100.0

    for spotFXRate in spotFXRates:
        callOption = FinFXVanillaOption(
            expiryDate, strikeFXRate, FinFXOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            spotFXRate,
            domDiscountCurve,
            forDiscountCurve,
            model)
        impliedVol = callOption.impliedVolatility(
            valueDate, spotFXRate, domDiscountCurve, forDiscountCurve, value)
        testCases.print(spotFXRate, value, volatility, impliedVol)


test_FinFXVanillaOption()
testCases.compareTestCases()
