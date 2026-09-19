# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path

_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause

_install_double_click_pause()
import add_fp_to_path

from financepy.products.equity.equity_one_touch_option import (
    EquityOneTouchOption,
)
from financepy.utils.global_types import TouchOptionTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

########################################################################################


def test_equity_one_touch_option():

    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    value_dt = Date(1, 1, 2016)
    expiry_dt = Date(2, 7, 2016)
    interest_rate = 0.10
    volatility = 0.20
    barrier_level = 100.0  # H
    model = BlackScholes(volatility)
    dividend_yield = 0.03
    num_paths = 10000
    num_steps_per_year = 252

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    stock_price = 105.0
    payment_size = 15.0

    print("================================= CASH ONLY")

    down_types = [
        TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT,
        TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
        TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING,
    ]

    print("TYPE", "VALUE", "VALUE_MC")

    for down_type in down_types:

        option = EquityOneTouchOption(expiry_dt, down_type, barrier_level, payment_size)

        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_mc = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_steps_per_year,
            num_paths,
        )

        print("%60s " % down_type, "%9.5f" % v, "%9.5f" % v_mc)

    stock_price = 95.0
    payment_size = 15.0

    up_types = [
        TouchOptionTypes.UP_AND_IN_CASH_AT_HIT,
        TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
        TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING,
    ]

    print("TYPE", "VALUE", "VALUE_MC")

    for up_type in up_types:

        option = EquityOneTouchOption(expiry_dt, up_type, barrier_level, payment_size)

        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_mc = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_steps_per_year,
            num_paths,
        )

        print("%60s " % up_type, "%9.5f" % v, "%9.5f" % v_mc)

    stock_price = 105.0

    print("================= ASSET ONLY")

    down_types = [
        TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT,
        TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
        TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING,
    ]

    print("TYPE", "VALUE", "VALUE_MC")
    for down_type in down_types:

        option = EquityOneTouchOption(expiry_dt, down_type, barrier_level)

        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_mc = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_steps_per_year,
            num_paths,
        )

        print("%60s " % down_type, "%9.5f" % v, "%9.5f" % v_mc)

    stock_price = 95.0

    up_types = [
        TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT,
        TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY,
        TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING,
    ]

    for up_type in up_types:

        option = EquityOneTouchOption(expiry_dt, up_type, barrier_level)

        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_mc = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_steps_per_year,
            num_paths,
        )

        print("%60s " % up_type, "%9.5f" % v, "%9.5f" % v_mc)


########################################################################################

test_equity_one_touch_option()
