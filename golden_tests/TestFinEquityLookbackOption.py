###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import time
from financepy.products.equity.equity_float_lookback_option import (
    EquityFloatLookbackOption,
)
from financepy.products.equity.equity_fixed_lookback_option import (
    EquityFixedLookbackOption,
)
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityLookBackOption():
    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.3
    interest_rate = 0.05
    dividend_yield = 0.01
    num_paths_range = [10000]
    stock_price_range = range(90, 110, 10)
    num_steps_per_year = 252

    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    ###############################################################################

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFloatLookbackOption(expiry_dt, option_type)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFloatLookbackOption(expiry_dt, option_type)
            stock_min = stock_price - 10
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFloatLookbackOption(expiry_dt, option_type)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFloatLookbackOption(expiry_dt, option_type)
            stock_max = stock_price + 10
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    ###############################################################################
    ###############################################################################

    stock_price_range = range(90, 110, 10)
    num_steps_per_year = 252

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_CALL
    k = 95.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_CALL
    k = 100.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_CALL
    k = 105.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price + 10.0
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_max,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_PUT
    k = 95.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_PUT
    k = 100.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    test_cases.header(
        "NUMPATHS",
        "OPTION_TYPE",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    option_type = OptionTypes.EUROPEAN_PUT
    k = 105.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = EquityFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price - 10.0
            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                volatility,
                stock_min,
                num_paths,
                num_steps_per_year,
            )
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value
            test_cases.print(
                num_paths,
                option_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )


###############################################################################


def test_example():

    expiry_dt = Date(1, 1, 2021)
    strike_price = 105.0
    option_typeCall = OptionTypes.EUROPEAN_CALL
    option_typePut = OptionTypes.EUROPEAN_PUT
    lookbackCall = EquityFixedLookbackOption(
        expiry_dt, option_typeCall, strike_price
    )
    lookbackPut = EquityFixedLookbackOption(
        expiry_dt, option_typePut, strike_price
    )

    value_dt = Date(1, 1, 2020)
    interest_rate = 0.10
    stock_price = 100.0
    dividend_yield = 0.0
    stock_min_max = 100.0

    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    volatilities = [0.30]

    test_cases.header("VALUE")
    for vol in volatilities:
        v = lookbackCall.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            vol,
            stock_min_max,
        )
        test_cases.print(v)


###############################################################################


test_example()
# test_EquityLookBackOption()
test_cases.compareTestCases()
