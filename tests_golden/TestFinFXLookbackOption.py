###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_float_lookback_option import FXFloatLookbackOption
from financepy.products.fx.fx_fixed_lookback_option import FXFixedLookbackOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityLookBackOption():
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.3
    num_pathsRange = [10000]
    stock_priceRange = range(90, 110, 5)
    num_steps_per_year = 252

    domesticRate = 0.05
    domestic_curve = DiscountCurveFlat(valuation_date, domesticRate)

    foreignRate = 0.02
    foreign_curve = DiscountCurveFlat(valuation_date, foreignRate)

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

    option_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFloatLookbackOption(expiry_date, option_type)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                stockMin,
                value,
                value_mc,
                diff,
                time_elapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFloatLookbackOption(expiry_date, option_type)
            stockMin = stock_price - 10
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                stockMin,
                value,
                value_mc,
                diff,
                time_elapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFloatLookbackOption(expiry_date, option_type)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                stockMax,
                value,
                value_mc,
                diff,
                time_elapsed)

    testCases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME")

    option_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFloatLookbackOption(expiry_date, option_type)
            stockMax = stock_price + 10
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                stockMax,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_CALL
    k = 95.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_CALL
    k = 100.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_CALL
    k = 105.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price + 10.0
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMax,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMax,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_PUT
    k = 95.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_PUT
    k = 100.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                time_elapsed)

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

    option_type = OptionTypes.EUROPEAN_PUT
    k = 105.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = FXFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price - 10.0
            value = option.value(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stockMin,
                num_paths,
                num_steps_per_year)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            testCases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stockMin,
                value,
                value_mc,
                diff,
                time_elapsed)

###############################################################################


test_EquityLookBackOption()
testCases.compareTestCases()
