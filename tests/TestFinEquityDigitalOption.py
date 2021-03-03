###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityDigitalOption import FinEquityDigitalOption, FinDigitalOptionTypes
from financepy.models.black_scholes import FinModelBlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityDigitalOption():

    underlyingType = FinDigitalOptionTypes.CASH_OR_NOTHING

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)
    
    model = FinModelBlackScholes(volatility)
    import time

    callOptionValues = []
    callOptionValuesMC = []
    num_pathsList = [
        10000,
        20000,
        40000,
        80000,
        160000,
        320000,
        640000,
        1280000,
        2560000]

    testCases.header("NumLoops", "ValueBS", "ValueMC", "TIME")

    for num_paths in num_pathsList:

        callOption = FinEquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, underlyingType)
        value = callOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        start = time.time()
        valueMC = callOption.valueMC(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model,
            num_paths)
        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, valueMC, duration)

        callOptionValues.append(value)
        callOptionValuesMC.append(valueMC)

#    plt.figure(figsize=(10,8))
#    plt.plot(num_pathsList, callOptionValues, color = 'b', label="Call Option")
#    plt.plot(num_pathsList, callOptionValuesMC, color = 'r', label = "Call Option MC")
#    plt.xlabel("Num Loops")
#    plt.legend(loc='best')

##########################################################################

    stock_prices = range(50, 150, 50)
    callOptionValues = []
    callOptionDeltas = []
    callOptionVegas = []
    callOptionThetas = []

    for stock_price in stock_prices:
        callOption = FinEquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, underlyingType)
        value = callOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        delta = callOption.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        vega = callOption.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        theta = callOption.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        callOptionValues.append(value)
        callOptionDeltas.append(delta)
        callOptionVegas.append(vega)
        callOptionThetas.append(theta)

    putOptionValues = []
    putOptionDeltas = []
    putOptionVegas = []
    putOptionThetas = []

    for stock_price in stock_prices:
        putOption = FinEquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, underlyingType)
        value = putOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        delta = putOption.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        vega = putOption.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        theta = putOption.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividendCurve,
            model)
        putOptionValues.append(value)
        putOptionDeltas.append(delta)
        putOptionVegas.append(vega)
        putOptionThetas.append(theta)

##########################################################################

test_FinEquityDigitalOption()
testCases.compareTestCases()
