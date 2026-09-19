# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import time

import add_fp_to_path

from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_float_lookback_option import FXFloatLookbackOption
from financepy.products.fx.fx_fixed_lookback_option import FXFixedLookbackOption
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date


########################################################################################


def test_equity_look_back_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.3
    num_paths_range = [10000]
    stock_price_range = range(90, 110, 5)
    num_steps_per_year = 252

    domestic_rate = 0.05
    domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)

    foreign_rate = 0.02
    foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFloatLookbackOption(expiry_dt, opt_type)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_CALL
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFloatLookbackOption(expiry_dt, opt_type)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFloatLookbackOption(expiry_dt, opt_type)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_PUT
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFloatLookbackOption(expiry_dt, opt_type)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    stock_price_range = range(90, 110, 5)
    num_steps_per_year = 252

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_CALL
    k = 95.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_CALL
    k = 100.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMAX",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_CALL
    k = 105.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_max,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_PUT
    k = 95.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_PUT
    k = 100.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

    print(
        "NUMPATHS",
        "opt_type",
        "S",
        "K",
        "SMIN",
        "VALUE",
        "VALUE_MC",
        "DIFF",
        "TIME",
    )

    opt_type = OptionTypes.EUROPEAN_PUT
    k = 105.0
    for stock_price in stock_price_range:
        for num_paths in num_paths_range:
            option = FXFixedLookbackOption(expiry_dt, opt_type, k)
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
            print(
                num_paths,
                opt_type,
                stock_price,
                k,
                stock_min,
                value,
                value_mc,
                diff,
                time_elapsed,
            )


########################################################################################

test_equity_look_back_option()
