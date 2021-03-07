###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.products.equity.equity_fixed_lookback_option import EquityFixedLookbackOption
from financepy.utils.global_types import FinOptionTypes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_EquityLookBackOption():
    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    stockPrice = 100.0
    volatility = 0.3
    interestRate = 0.05
    dividendYield = 0.01
    numPathsRange = [10000]
    stockPriceRange = range(90, 110, 10)
    numStepsPerYear = 252

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

###############################################################################

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_CALL
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFloatLookbackOption(expiryDate, option_type)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                stockMin,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_CALL
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFloatLookbackOption(expiryDate, option_type)
            stockMin = stockPrice - 10
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                stockMin,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_PUT
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFloatLookbackOption(expiryDate, option_type)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                stockMax,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_PUT
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFloatLookbackOption(expiryDate, option_type)
            stockMax = stockPrice + 10
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                stockMax,
                value,
                value_mc,
                diff,
                timeElapsed)

###############################################################################
###############################################################################

    stockPriceRange = range(90, 110, 10)
    numStepsPerYear = 252

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_CALL
    k = 95.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_CALL
    k = 100.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_CALL
    k = 105.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice + 10.0
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMax,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_PUT
    k = 95.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_PUT
    k = 100.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                timeElapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = FinOptionTypes.EUROPEAN_PUT
    k = 105.0
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = EquityFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice - 10.0
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                volatility,
                stockMin,
                numPaths,
                numStepsPerYear)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                numPaths,
                option_type,
                stockPrice,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                timeElapsed)

###############################################################################

def test_example():

    expiryDate = FinDate(1, 1, 2021)
    strike_price = 105.0
    option_typeCall = FinOptionTypes.EUROPEAN_CALL
    option_typePut = FinOptionTypes.EUROPEAN_PUT
    lookbackCall = EquityFixedLookbackOption(expiryDate, option_typeCall, strike_price)
    lookbackPut = EquityFixedLookbackOption(expiryDate, option_typePut, strike_price)

    valueDate = FinDate(1, 1, 2020)
    interestRate = 0.10
    stockPrice = 100.0
    dividendYield = 0.0
    stockMinMax = 100.0

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    volatilities = [0.30]

    testCases.header("VALUE")
    for vol in volatilities:
        v = lookbackCall.value(valueDate, stockPrice, discountCurve, dividend_curve, vol, stockMinMax)
        testCases.print(v)

###############################################################################

test_example()
#test_EquityLookBackOption()
testCases.compareTestCases()
