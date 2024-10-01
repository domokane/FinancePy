###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import time
from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_float_lookback_option import (
    FXFloatLookbackOption,
)
from financepy.products.fx.fx_fixed_lookback_option import (
    FXFixedLookbackOption,
)
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
    num_paths_range = [10000]
    stock_price_range = range(90, 110, 5)
    num_steps_per_year = 252

    domesticRate = 0.05
    domestic_curve = DiscountCurveFlat(value_dt, domesticRate)

    foreignRate = 0.02
    foreign_curve = DiscountCurveFlat(value_dt, foreignRate)

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
            option = FXFloatLookbackOption(expiry_dt, option_type)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFloatLookbackOption(expiry_dt, option_type)
            stock_min = stock_price - 10
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFloatLookbackOption(expiry_dt, option_type)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFloatLookbackOption(expiry_dt, option_type)
            stock_max = stock_price + 10
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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

    stock_price_range = range(90, 110, 5)
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_max = stock_price + 10.0
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_max,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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
            option = FXFixedLookbackOption(expiry_dt, option_type, k)
            stock_min = stock_price - 10.0
            value = option.value(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
                volatility,
                stock_min,
            )
            start = time.time()
            value_mc = option.value_mc(
                value_dt,
                stock_price,
                domestic_curve,
                foreign_curve,
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


test_EquityLookBackOption()
test_cases.compareTestCases()
