# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import time
import numpy as np

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.FinOptionTypes import FinOptionTypes
from financepy.products.fx.FinFXForward import FinFXForward
from financepy.products.fx.FinFXVanillaOption import FinFXVanillaOption
from financepy.products.fx.FinFXModelTypes import FinFXModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve

from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.market.curves.FinLiborCurve import FinLiborCurve
from financepy.products.libor.FinLiborDeposit import FinLiborDeposit

from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################

def test_FinFXVanillaOptionWystupExample1():

    # Example from Book extract by Uwe Wystup with results in Table 1.2
    # https://mathfinance.com/wp-content/uploads/2017/06/FXOptionsStructuredProducts2e-Extract.pdf

    # Not exactly T=1.0 but close so don't exact exact agreement
    # (in fact I do not get exact agreement even if I do set T=1.0)
    valueDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.030  # EUR
    ccy2CCRate = 0.025  # USD

    currencyPair = ccy1 + ccy2  # Always ccy1ccy2
    spotFXRate = 1.20
    strikeFXRate = 1.250
    volatility = 0.10

    spotDays = 0
    settlementDate = valueDate.addWorkDays(spotDays)
    maturityDate = settlementDate.addMonths(12)
    notional = 1000000.0
    notionalCurrency = "EUR"
    calendarType = FinCalendarTypes.TARGET

    domDiscountCurve = FinFlatCurve(valueDate, ccy2CCRate)
    forDiscountCurve = FinFlatCurve(valueDate, ccy1CCRate)

    model = FinFXModelBlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    callOption = FinFXVanillaOption(expiryDate,
                                    strikeFXRate,
                                    currencyPair,
                                    FinOptionTypes.EUROPEAN_CALL,
                                    notional,
                                    "EUR", 2)

    value = callOption.value(
        1.0,
        spotFXRate,
        domDiscountCurve,
        forDiscountCurve,
        model)

    notional = 1250000.0
    callOption = FinFXVanillaOption(expiryDate,
                                    strikeFXRate,
                                    currencyPair,
                                    FinOptionTypes.EUROPEAN_CALL,
                                    notional,
                                    "USD", 2)

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

    testCases.header("value", "delta")
    testCases.print(value, delta)

###############################################################################

def test_FinFXVanillaOptionWystupExample2():

    # Example Bloomberg Pricing at
    # https://stackoverflow.com/questions/48778712/fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valueDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.0396  # EUR
    ccy2CCRate = 0.0357  # USD

    currencyPair = ccy1 + ccy2  # Always ccy1ccy2
    spotFXRate = 0.9090
    strikeFXRate = 0.9090
    volatility = 0.12

    spotDays = 0
    settlementDate = valueDate.addWorkDays(spotDays)
    maturityDate = settlementDate.addMonths(12)
    notional = 1000000.0
    notionalCurrency = "EUR"
    calendarType = FinCalendarTypes.TARGET

    domDiscountCurve = FinFlatCurve(valueDate, ccy2CCRate)
    forDiscountCurve = FinFlatCurve(valueDate, ccy1CCRate)

    model = FinFXModelBlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    callOption = FinFXVanillaOption(expiryDate,
                                    strikeFXRate,
                                    currencyPair,
                                    FinOptionTypes.EUROPEAN_PUT,
                                    notional,
                                    "EUR", 2)

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

    testCases.header("value", "delta")
    testCases.print(value, delta)

###############################################################################

def test_FinFXVanillaOptionBloombergExample():

    # Example Bloomberg Pricing at
    # https://stackoverflow.com/questions/48778712/fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valueDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(15, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    forName = "EUR"
    domName = "USD"
    forDepoRate = 0.05  # EUR
    domDepoRate = 0.02  # USD

    currencyPair = forName + domName  # Always FORDOM
    spotFXRate = 1.30
    strikeFXRate = 1.3650
    volatility = 0.20

    spotDays = 0
    settlementDate = valueDate.addWorkDays(spotDays)
    maturityDate = settlementDate.addMonths(12)
    notional = 1000000.0
    notionalCurrency = "EUR"
    calendarType = FinCalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    depo = FinLiborDeposit(settlementDate, maturityDate, domDepoRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    domDiscountCurve = FinLiborCurve(forName, settlementDate, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    depo = FinLiborDeposit(settlementDate, maturityDate, forDepoRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    forDiscountCurve = FinLiborCurve(domName, settlementDate, depos, fras, swaps)

    model = FinFXModelBlackScholes(volatility)

    callOption = FinFXVanillaOption(expiryDate,
                                    strikeFXRate,
                                    currencyPair,
                                    FinOptionTypes.EUROPEAN_CALL,
                                    notional,
                                    notionalCurrency, 2)

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

    testCases.header("value", "delta")
    testCases.print(value, delta)

###############################################################################

def test_FinFXVanillaOptionHullExample():

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

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_CALL,
                                        1000000,
                                        "USD")

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

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_CALL,
                                        1000000,
                                        "USD")

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

        putOption = FinFXVanillaOption(expiryDate,
                                       strikeFXRate,
                                       "EURUSD",
                                       FinOptionTypes.EUROPEAN_PUT,
                                       1000000,
                                       "USD")

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
        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_CALL,
                                        1000000,
                                        "USD")
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
        putOption = FinFXVanillaOption(expiryDate,
                                       strikeFXRate,
                                       "EURUSD",
                                       FinOptionTypes.EUROPEAN_PUT,
                                       1000000,
                                       "USD")

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
        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_CALL,
                                        1000000,
                                        "USD")

        value = callOption.value(valueDate,
                                 spotFXRate,
                                 domDiscountCurve,
                                 forDiscountCurve,
                                 model)['v']

        impliedVol = callOption.impliedVolatility(valueDate,
                                                  spotFXRate,
                                                  domDiscountCurve,
                                                  forDiscountCurve,
                                                  value)

        testCases.print(spotFXRate, value, volatility, impliedVol)

###############################################################################

test_FinFXVanillaOptionWystupExample1()
test_FinFXVanillaOptionWystupExample2()
test_FinFXVanillaOptionBloombergExample()
test_FinFXVanillaOptionHullExample()
testCases.compareTestCases()
