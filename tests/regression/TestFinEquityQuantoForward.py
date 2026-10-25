# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import math
import add_fp_to_path

from financepy.products.equity.equity_quanto_forward import EquityQuantoForward
from financepy.utils.date import Date
from financepy.utils.global_types import LongShortTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_quanto_forward():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    forward_price = 125.0  # Locked strike

    domestic_discount_rate = 0.05
    foreign_discount_rate = 0.03
    dividend_rate = 0.02

    equity_volatility = 0.25
    fx_volatility = 0.12
    correlation = 0.50

    notional = 100.0

    domestic_discount_curve = FlatDiscountCurve(
        value_dt, domestic_discount_rate
    )

    foreign_discount_curve = FlatDiscountCurve(
        value_dt, foreign_discount_rate
    )

    dividend_curve = FlatDiscountCurve(
        value_dt, dividend_rate
    )

    equity_quanto_forward = EquityQuantoForward(
        expiry_dt,
        forward_price,
        notional,
        LongShortTypes.LONG,
    )

    test_cases.header(
        "STOCK PRICE",
        "QUANTO FWD",
        "QUANTO VALUE",
        "QUANTO ADJ",
    )

    fwd_price = equity_quanto_forward.forward(
        value_dt,
        stock_price,
        foreign_discount_curve,
        dividend_curve,
        equity_volatility,
        fx_volatility,
        correlation,
    )

    fwd_value = equity_quanto_forward.value(
        value_dt,
        stock_price,
        domestic_discount_curve,
        foreign_discount_curve,
        dividend_curve,
        equity_volatility,
        fx_volatility,
        correlation,
    )

    tau = (expiry_dt - value_dt) / 365.0

    quanto_adjustment = math.exp(
        -correlation * equity_volatility * fx_volatility * tau
    )

    test_cases.print(
        stock_price,
        fwd_price,
        fwd_value,
        quanto_adjustment,
    )


########################################################################################


def test_zero_correlation():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    forward_price = 125.0
    notional = 100.0

    foreign_curve = FlatDiscountCurve(value_dt, 0.03)
    dividend_curve = FlatDiscountCurve(value_dt, 0.02)

    equity_quanto_forward = EquityQuantoForward(
        expiry_dt,
        forward_price,
        notional,
        LongShortTypes.LONG,
    )

    fwd_price = equity_quanto_forward.forward(
        value_dt,
        stock_price,
        foreign_curve,
        dividend_curve,
        0.25,
        0.12,
        0.0,
    )

    expected_forward = (
        stock_price
        * dividend_curve.df(expiry_dt)
        / foreign_curve.df(expiry_dt)
    )

    assert math.isclose(
        fwd_price, expected_forward, rel_tol=1e-10
    )

    test_cases.header("ZERO CORR FWD", "EXPECTED FWD")
    test_cases.print(fwd_price, expected_forward)


########################################################################################


def test_long_short_symmetry():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    domestic_curve = FlatDiscountCurve(value_dt, 0.05)
    foreign_curve = FlatDiscountCurve(value_dt, 0.03)
    dividend_curve = FlatDiscountCurve(value_dt, 0.02)

    long_forward = EquityQuantoForward(
        expiry_dt, 125.0, 100.0, LongShortTypes.LONG
    )

    short_forward = EquityQuantoForward(
        expiry_dt, 125.0, 100.0, LongShortTypes.SHORT
    )

    args = (
        value_dt,
        130.0,
        domestic_curve,
        foreign_curve,
        dividend_curve,
        0.25,
        0.12,
        0.50,
    )

    long_value = long_forward.value(*args)
    short_value = short_forward.value(*args)

    assert math.isclose(
        long_value, -short_value, rel_tol=1e-10
    )

    test_cases.header("LONG VALUE", "SHORT VALUE")
    test_cases.print(long_value, short_value)


########################################################################################

test_equity_quanto_forward()
test_zero_correlation()
test_long_short_symmetry()

test_cases.compare_test_cases()
