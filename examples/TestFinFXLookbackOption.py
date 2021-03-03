###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.products.fx.FinFXFloatLookbackOption import FinFXFloatLookbackOption
from financepy.products.fx.FinFXFixedLookbackOption import FinFXFixedLookbackOption
from financepy.market.curves.FinDiscountCurveFlat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinEquityLookBackOption():
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.3
    num_pathsRange = [10000]
    stock_priceRange = range(90, 110, 5)
    num_steps_per_year = 252

    domesticRate = 0.05
    domesticCurve = DiscountCurveFlat(valuation_date, domesticRate)

    foreignRate = 0.02
    foreignCurve = DiscountCurveFlat(valuation_date, foreignRate)

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
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMin = stock_price - 10
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMax = stock_price + 10
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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

    stock_priceRange = range(90, 110, 5)
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stock_price + 10.0
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stock_price - 10.0
            value = option.value(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                domesticCurve,
                foreignCurve,
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


test_FinEquityLookBackOption()
testCases.compareTestCases()
