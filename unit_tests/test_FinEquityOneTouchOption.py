# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_types import TouchOptionTypes
from financepy.products.equity.equity_one_touch_option import (
    EquityOneTouchOption,
)

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

payment_size = 15.0


# Made error bar larger as values larger
def assert_close(value, expected, tol=2.5e-2):
    assert np.isclose(value, expected, atol=tol), (
        f"value={value:.10f}, expected={expected:.10f}, "
        f"diff={value - expected:.10f}"
    )


########################################################################################


def test_down_and_in_cash_at_hit():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT
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

    assert_close(v, 10.15381)
    assert_close(v_mc, 9.88882)


########################################################################################


def test_down_and_in_cash_at_expiry():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY
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

    assert_close(v, 9.77218)
    assert_close(v_mc, 9.51229)


########################################################################################


def test_down_and_out_cash_or_nothing():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING
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

    assert_close(v, 4.49627)
    assert_close(v_mc, 4.75615)


########################################################################################


def test_up_and_in_cash_at_hit():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_IN_CASH_AT_HIT
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

    assert_close(v, 11.28531)
    assert_close(v_mc, 11.11317)


########################################################################################


def test_up_and_in_cash_at_expiry():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY
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

    assert_close(v, 10.86668)
    assert_close(v_mc, 10.70133)


########################################################################################


def test_up_and_out_cash_or_nothing():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING
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

    assert_close(v, 3.40176)
    assert_close(v_mc, 3.56711)


########################################################################################


def test_down_and_in_asset_at_hit():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT
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

    assert_close(v, 67.69205)
    assert_close(v_mc, 65.92547)


########################################################################################


def test_down_and_in_asset_at_expiry():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY
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

    assert_close(v, 66.91760)
    assert_close(v_mc, 66.66667)


########################################################################################


def test_down_and_out_asset_or_nothing():

    stock_price = 105.0
    down_type = TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING
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

    assert_close(v, 36.51916)
    assert_close(v_mc, 38.66099)


########################################################################################


def test_up_and_in_asset_at_hit():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT
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

    assert_close(v, 75.23538)
    assert_close(v_mc, 74.08783)


########################################################################################


def test_up_and_in_asset_at_expiry():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY
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

    assert_close(v, 74.38596)
    assert_close(v_mc, 75.00000)


########################################################################################


def test_up_and_out_asset_or_nothing():

    stock_price = 95.0
    down_type = TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING
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

    assert_close(v, 19.19968)
    assert_close(v_mc, 20.00701)
