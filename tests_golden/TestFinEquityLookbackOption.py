###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityFloatLookbackOption import FinEquityFloatLookbackOption
from financepy.products.equity.FinEquityFixedLookbackOption import FinEquityFixedLookbackOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinEquityLookBackOption():
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.3
    interestRate = 0.05
    dividendYield = 0.01
    num_pathsRange = [10000]
    stock_priceRange = range(90, 110, 10)
    num_steps_per_year = 252

    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

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

    optionType = FinOptionTypes.EUROPEAN_CALL
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFloatLookbackOption(expiry_date, optionType)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                stockMin,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_CALL
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFloatLookbackOption(expiry_date, optionType)
            stockMin = stock_price - 10
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                stockMin,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_PUT
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFloatLookbackOption(expiry_date, optionType)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                stockMax,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_PUT
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFloatLookbackOption(expiry_date, optionType)
            stockMax = stock_price + 10
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                stockMax,
                value,
                valueMC,
                diff,
                timeElapsed)

###############################################################################
###############################################################################

    stock_priceRange = range(90, 110, 10)
    num_steps_per_year = 252

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

    optionType = FinOptionTypes.EUROPEAN_CALL
    k = 95.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMax,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_CALL
    k = 100.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMax,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_CALL
    k = 105.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price + 10.0
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMax,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_PUT
    k = 95.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMin,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_PUT
    k = 100.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMin,
                value,
                valueMC,
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

    optionType = FinOptionTypes.EUROPEAN_PUT
    k = 105.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FinEquityFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price - 10.0
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                num_paths,
                optionType,
                stock_price,
                k,
                stockMin,
                value,
                valueMC,
                diff,
                timeElapsed)

###############################################################################

def test_example():

    expiry_date = Date(1, 1, 2021)
    strikePrice = 105.0
    optionTypeCall = FinOptionTypes.EUROPEAN_CALL
    optionTypePut = FinOptionTypes.EUROPEAN_PUT
    lookbackCall = FinEquityFixedLookbackOption(expiry_date, optionTypeCall, strikePrice)
    lookbackPut = FinEquityFixedLookbackOption(expiry_date, optionTypePut, strikePrice)

    valuation_date = Date(1, 1, 2020)
    interestRate = 0.10
    stock_price = 100.0
    dividendYield = 0.0
    stockMinMax = 100.0

    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    volatilities = [0.30]

    testCases.header("VALUE")
    for vol in volatilities:
        v = lookbackCall.value(valuation_date, stock_price, discount_curve, dividendCurve, vol, stockMinMax)
        testCases.print(v)

###############################################################################

test_example()
#test_FinEquityLookBackOption()
testCases.compareTestCases()
