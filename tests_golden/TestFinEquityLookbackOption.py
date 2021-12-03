###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.products.equity.equity_fixed_lookback_option import EquityFixedLookbackOption
from financepy.utils.global_types import OptionTypes
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
    interest_rate = 0.05
    dividend_yield = 0.01
    num_pathsRange = [10000]
    stock_priceRange = range(90, 110, 10)
    num_steps_per_year = 252

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

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
            option = EquityFloatLookbackOption(expiry_date, option_type)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFloatLookbackOption(expiry_date, option_type)
            stockMin = stock_price - 10
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFloatLookbackOption(expiry_date, option_type)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFloatLookbackOption(expiry_date, option_type)
            stockMax = stock_price + 10
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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

    option_type = OptionTypes.EUROPEAN_CALL
    k = 95.0
    for stock_price in stock_priceRange:
        for num_paths in num_pathsRange:
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMax = stock_price + 10.0
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMax)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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
            option = EquityFixedLookbackOption(expiry_date, option_type, k)
            stockMin = stock_price - 10.0
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stockMin)
            start = time.time()
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
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


def test_example():

    expiry_date = Date(1, 1, 2021)
    strike_price = 105.0
    option_typeCall = OptionTypes.EUROPEAN_CALL
    option_typePut = OptionTypes.EUROPEAN_PUT
    lookbackCall = EquityFixedLookbackOption(
        expiry_date, option_typeCall, strike_price)
    lookbackPut = EquityFixedLookbackOption(
        expiry_date, option_typePut, strike_price)

    valuation_date = Date(1, 1, 2020)
    interest_rate = 0.10
    stock_price = 100.0
    dividend_yield = 0.0
    stock_min_max = 100.0

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    volatilities = [0.30]

    testCases.header("VALUE")
    for vol in volatilities:
        v = lookbackCall.value(valuation_date, stock_price,
                               discount_curve, dividend_curve, vol, stock_min_max)
        testCases.print(v)

###############################################################################


test_example()
# test_EquityLookBackOption()
testCases.compareTestCases()
