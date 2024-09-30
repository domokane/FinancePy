###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import TouchOptionTypes
from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption

value_dt = Date(1, 1, 2016)
expiry_dt = Date(2, 7, 2016)
interest_rate = 0.10
volatility = 0.20
barrier_level = 100.0  # H
model = BlackScholes(volatility)
dividend_yield = 0.03
num_paths = 10000
num_steps_per_year = 252

discount_curve = DiscountCurveFlat(value_dt, interest_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

payment_size = 15.0


def test_DOWN_AND_IN_CASH_AT_HIT():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 10.15381
    assert round(v_mc, 5) == 9.88882


def test_DOWN_AND_IN_CASH_AT_EXPIRY():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 9.77218
    assert round(v_mc, 5) == 9.51229


def test_DOWN_AND_OUT_CASH_OR_NOTHING():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 4.49627
    assert round(v_mc, 5) == 4.75615


def test_UP_AND_IN_CASH_AT_HIT():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_IN_CASH_AT_HIT
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 11.28531
    assert round(v_mc, 5) == 11.11317


def test_UP_AND_IN_CASH_AT_EXPIRY():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 10.86668
    assert round(v_mc, 5) == 10.70133


def test_UP_AND_OUT_CASH_OR_NOTHING():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 3.40176
    assert round(v_mc, 5) == 3.56711


def test_DOWN_AND_IN_ASSET_AT_HIT():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 67.69205
    assert round(v_mc, 5) == 65.92547


def test_DOWN_AND_IN_ASSET_AT_EXPIRY():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 66.91760
    assert round(v_mc, 5) == 66.66667


def test_DOWN_AND_OUT_ASSET_OR_NOTHING():
    stock_price = 105.0
    downType = TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 36.51916
    assert round(v_mc, 5) == 38.66099


def test_UP_AND_IN_ASSET_AT_HIT():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 75.23538
    assert round(v_mc, 5) == 74.08783


def test_UP_AND_IN_ASSET_AT_EXPIRY():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 74.38596
    assert round(v_mc, 5) == 75.00000


def test_UP_AND_OUT_ASSET_OR_NOTHING():
    stock_price = 95.0
    downType = TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING
    option = EquityOneTouchOption(expiry_dt, downType, barrier_level, payment_size)
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

    assert round(v, 5) == 19.19968
    assert round(v_mc, 5) == 20.00701
