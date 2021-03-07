###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
import sys
sys.path.append("..")

from financepy.utils.global_types import FinOptionTypes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.FinIborSingleCurve import IborSingleCurve
from financepy.products.rates.FinIborDeposit import FinIborDeposit
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
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

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 1.20
    strike_fx_rate = 1.250
    volatility = 0.10

    notional = 1000000.0

    dom_discount_curve = FinDiscountCurveFlat(valueDate, ccy2CCRate)
    for_discount_curve = FinDiscountCurveFlat(valueDate, ccy1CCRate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    callOption = FXVanillaOption(expiryDate,
                                 strike_fx_rate,
                                 currency_pair,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 notional,
                                    "EUR", 2)

    value = callOption.value(
        1.0,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    notional = 1250000.0
    callOption = FXVanillaOption(expiryDate,
                                 strike_fx_rate,
                                 currency_pair,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 notional,
                                    "USD", 2)

    value = callOption.value(
        valueDate,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    delta = callOption.delta(
        valueDate,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
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

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 0.9090
    strike_fx_rate = 0.9090
    volatility = 0.12

    notional = 1000000.0

    dom_discount_curve = FinDiscountCurveFlat(valueDate, ccy2CCRate)
    for_discount_curve = FinDiscountCurveFlat(valueDate, ccy1CCRate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0
    callOption = FXVanillaOption(expiryDate,
                                 strike_fx_rate,
                                 currency_pair,
                                 FinOptionTypes.EUROPEAN_PUT,
                                 notional,
                                    "EUR", 2)

    value = callOption.value(
            valueDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)

    delta = callOption.delta(
            valueDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)

    testCases.header("value", "delta")
    testCases.print(value, delta)

###############################################################################

def test_FinFXVanillaOptionBloombergExample():

    # Example Bloomberg Pricing at
    # https://stackoverflow.com/questions/48778712/fx-vanilla-call-price-in-quantlib-doesnt-match-bloomberg

    valuationDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(15, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    forName = "EUR"
    domName = "USD"
    forDepoRate = 0.05  # EUR
    domDepoRate = 0.02  # USD

    currency_pair = forName + domName  # Always FORDOM
    spot_fx_rate = 1.30
    strike_fx_rate = 1.3650
    volatility = 0.20

    spot_days = 0
    settlementDate = valuationDate.addWeekDays(spot_days)
    maturityDate = settlementDate.addMonths(12)
    notional = 1000000.0
    notional_currency = "EUR"
    calendarType = FinCalendarTypes.TARGET

    depos = []
    fras = []
    swaps = []
    depo = FinIborDeposit(settlementDate, maturityDate, domDepoRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    dom_discount_curve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    depos = []
    fras = []
    swaps = []
    depo = FinIborDeposit(settlementDate, maturityDate, forDepoRate,
                           FinDayCountTypes.ACT_360, notional, calendarType)
    depos.append(depo)
    for_discount_curve = FinIborSingleCurve(valuationDate, depos, fras, swaps)

    model = BlackScholes(volatility)

    callOption = FXVanillaOption(expiryDate,
                                 strike_fx_rate,
                                 currency_pair,
                                 FinOptionTypes.EUROPEAN_CALL,
                                 notional,
                                 notional_currency, 2)

    value = callOption.value(
        valuationDate,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    delta = callOption.delta(
        valuationDate,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

    testCases.header("value", "delta")
    testCases.print(value, delta)

###############################################################################


def test_FinFXVanillaOptionHullExample():

    #   Example from Hull 4th edition page 284
    valuationDate = FinDate(1, 1, 2015)
    expiryDate = valuationDate.addMonths(4)
    spot_fx_rate = 1.60
    volatility = 0.1411
    dom_interest_rate = 0.08
    forInterestRate = 0.11
    model = BlackScholes(volatility)
    dom_discount_curve = FinDiscountCurveFlat(valuationDate, dom_interest_rate)
    for_discount_curve = FinDiscountCurveFlat(valuationDate, forInterestRate)

    numPathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")
    strike_fx_rate = 1.60

    for numPaths in numPathsList:

        callOption = FXVanillaOption(expiryDate,
                                     strike_fx_rate,
                                        "EURUSD",
                                     FinOptionTypes.EUROPEAN_CALL,
                                     1000000,
                                        "USD")

        value = callOption.value(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)

        start = time.time()

        value_mc = callOption.value_mc(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model,
            numPaths)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, value_mc, duration)

##########################################################################

    spot_fx_rates = np.arange(100, 200, 10)
    spot_fx_rates = spot_fx_rates/100.0
    numPaths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC", "TIME")

    for spot_fx_rate in spot_fx_rates:

        callOption = FXVanillaOption(expiryDate,
                                     strike_fx_rate,
                                        "EURUSD",
                                     FinOptionTypes.EUROPEAN_CALL,
                                     1000000,
                                        "USD")

        value = callOption.value(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        start = time.time()
        value_mc = callOption.value_mc(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, value_mc, duration)

##########################################################################

    spot_fx_rates = np.arange(100, 200, 10) / 100.0
    numPaths = 100000

    testCases.header("SPOT FX RATE", "PUT_VALUE_BS", "PUT_VALUE_MC", "TIME")

    for spot_fx_rate in spot_fx_rates:

        putOption = FXVanillaOption(expiryDate,
                                    strike_fx_rate,
                                       "EURUSD",
                                    FinOptionTypes.EUROPEAN_PUT,
                                    1000000,
                                       "USD")

        value = putOption.value(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        start = time.time()
        value_mc = putOption.value_mc(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(spot_fx_rate, value, value_mc, duration)

##########################################################################

    spot_fx_rates = np.arange(100, 200, 10)/100.0

    testCases.header(
        "SPOT FX RATE",
        "CALL_VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for spot_fx_rate in spot_fx_rates:
        callOption = FXVanillaOption(expiryDate,
                                     strike_fx_rate,
                                        "EURUSD",
                                     FinOptionTypes.EUROPEAN_CALL,
                                     1000000,
                                        "USD")
        value = callOption.value(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        delta = callOption.delta(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        vega = callOption.vega(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        theta = callOption.theta(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        #  callOption.rho(valueDate,stockPrice, interestRate,
        #  dividendYield, modelType, model_params)
        rho = 999
        testCases.print(spot_fx_rate, value, delta, vega, theta, rho)

    testCases.header(
        "SPOT FX RATE",
        "PUT_VALUE_BS",
        "DELTA_BS",
        "VEGA_BS",
        "THETA_BS",
        "RHO_BS")

    for spot_fx_rate in spot_fx_rates:
        putOption = FXVanillaOption(expiryDate,
                                    strike_fx_rate,
                                       "EURUSD",
                                    FinOptionTypes.EUROPEAN_PUT,
                                    1000000,
                                       "USD")

        value = putOption.value(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        delta = putOption.delta(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        vega = putOption.vega(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        theta = putOption.theta(
            valuationDate,
            spot_fx_rate,
            dom_discount_curve,
            for_discount_curve,
            model)
        # putOption.rho(valueDate,stockPrice, interestRate, dividendYield,
        # modelType, model_params)
        rho = 999
        testCases.print(spot_fx_rate, value, delta, vega, theta, rho)

##########################################################################

    testCases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spot_fx_rates = np.arange(100, 200, 10)/100.0

    for spot_fx_rate in spot_fx_rates:
        callOption = FXVanillaOption(expiryDate,
                                     strike_fx_rate,
                                        "EURUSD",
                                     FinOptionTypes.EUROPEAN_CALL,
                                     1000000,
                                        "USD")

        value = callOption.value(valuationDate,
                                 spot_fx_rate,
                                 dom_discount_curve,
                                 for_discount_curve,
                                 model)['v']

        impliedVol = callOption.implied_volatility(valuationDate,
                                                  spot_fx_rate,
                                                  dom_discount_curve,
                                                  for_discount_curve,
                                                  value)

        testCases.print(spot_fx_rate, value, volatility, impliedVol)

###############################################################################


test_FinFXVanillaOptionWystupExample1()
test_FinFXVanillaOptionWystupExample2()
test_FinFXVanillaOptionBloombergExample()
test_FinFXVanillaOptionHullExample()
testCases.compareTestCases()
