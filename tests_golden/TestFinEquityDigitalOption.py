###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.market.curves.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_digital_option import EquityDigitalOption, FinDigitalOptionTypes
from financepy.utils.global_types import FinOptionTypes
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityDigitalOption():

    underlying_type = FinDigitalOptionTypes.CASH_OR_NOTHING

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    model = BlackScholes(volatility)
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

        callOption = EquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, underlying_type)
        value = callOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        start = time.time()
        value_mc = callOption.value_mc(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths)
        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc, duration)

        callOptionValues.append(value)
        callOptionValuesMC.append(value_mc)

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
        callOption = EquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL, underlying_type)
        value = callOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        delta = callOption.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        vega = callOption.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        theta = callOption.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
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
        putOption = EquityDigitalOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_PUT, underlying_type)
        value = putOption.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        delta = putOption.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        vega = putOption.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        theta = putOption.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        putOptionValues.append(value)
        putOptionDeltas.append(delta)
        putOptionVegas.append(vega)
        putOptionThetas.append(theta)

##########################################################################


test_EquityDigitalOption()
testCases.compareTestCases()
