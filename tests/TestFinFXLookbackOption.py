###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.utils.global_types import FinOptionTypes
from financepy.products.fx.fx_float_lookback_option import FinFXFloatLookbackOption
from financepy.products.fx.fx_fixed_lookback_option import FinFXFixedLookbackOption
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
    numPathsRange = [10000]
    stockPriceRange = range(90, 110, 5)
    numStepsPerYear = 252

    domesticRate = 0.05
    domestic_curve = FinDiscountCurveFlat(valueDate, domesticRate)

    foreignRate = 0.02
    foreign_curve = FinDiscountCurveFlat(valueDate, foreignRate)

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
            option = FinFXFloatLookbackOption(expiryDate, option_type)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFloatLookbackOption(expiryDate, option_type)
            stockMin = stockPrice - 10
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFloatLookbackOption(expiryDate, option_type)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFloatLookbackOption(expiryDate, option_type)
            stockMax = stockPrice + 10
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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

    stockPriceRange = range(90, 110, 5)
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMax = stockPrice + 10.0
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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
            option = FinFXFixedLookbackOption(expiryDate, option_type, k)
            stockMin = stockPrice - 10.0
            value = option.value(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                domestic_curve,
                foreign_curve,
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


test_EquityLookBackOption()
testCases.compareTestCases()
