###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.products.fx.FinFXFloatLookbackOption import FinFXFloatLookbackOption
from financepy.products.fx.FinFXFixedLookbackOption import FinFXFixedLookbackOption
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.utils.Date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinEquityLookBackOption():
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stockPrice = 100.0
    volatility = 0.3
    numPathsRange = [10000]
    stockPriceRange = range(90, 110, 5)
    num_steps_per_year = 252

    domesticRate = 0.05
    domesticCurve = FinDiscountCurveFlat(valuation_date, domesticRate)

    foreignRate = 0.02
    foreignCurve = FinDiscountCurveFlat(valuation_date, foreignRate)

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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMin = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMin = stockPrice - 10
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMax = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFloatLookbackOption(expiry_date, optionType)
            stockMax = stockPrice + 10
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
                stockMax,
                value,
                valueMC,
                diff,
                timeElapsed)

###############################################################################
###############################################################################

    stockPriceRange = range(90, 110, 5)
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMax = stockPrice + 10.0
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMax,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stockPrice
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
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
    for stockPrice in stockPriceRange:
        for numPaths in numPathsRange:
            option = FinFXFixedLookbackOption(expiry_date, optionType, k)
            stockMin = stockPrice - 10.0
            value = option.value(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin)
            start = time.time()
            valueMC = option.valueMC(
                valuation_date,
                stockPrice,
                domesticCurve,
                foreignCurve,
                volatility,
                stockMin,
                numPaths,
                num_steps_per_year)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value
            testCases.print(
                numPaths,
                optionType,
                stockPrice,
                k,
                stockMin,
                value,
                valueMC,
                diff,
                timeElapsed)

###############################################################################


test_FinEquityLookBackOption()
testCases.compareTestCases()
