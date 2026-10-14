########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
import pytest
from financepy.utils.error import FinError

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_european_analytic import european_value
from financepy.products.equity.equity_forward_start_option import (
    EquityForwardStartOption,
)

value_dt = Date(1, 1, 2026)
start_dt = Date(1, 7, 2026)
expiry_dt = Date(1, 7, 2027)
stock_price = 100.0
volatility = 0.25
interest_rate = 0.05
dividend_yield = 0.02
model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)


@pytest.mark.parametrize("spot", [0.0, -1.0, np.nan, np.inf, -np.inf, np.array([100.0, np.nan]), np.array([])])
def test_forward_start_rejects_invalid_stock_prices(spot):
    """Invalid stock levels fail with FinError before pricing cashflows."""
    option = EquityForwardStartOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL)
    with pytest.raises(FinError, match="Stock price"):
        option.value(value_dt, spot, discount_curve, dividend_curve, model)


def test_forward_start_preserves_valid_vector_price_scaling():
    """The validator preserves vector prices and homogeneous payoff scaling."""
    option = EquityForwardStartOption(start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL)
    scalar = option.value(value_dt, 100.0, discount_curve, dividend_curve, model)
    vector = option.value(value_dt, np.array([100.0, 200.0]), discount_curve, dividend_curve, model)
    np.testing.assert_allclose(vector, [scalar, 2 * scalar])


########################################################################################


def _rubinstein(opt_type, t_start, t_expiry):
    """Rubinstein (1991) at-the-money forward-start option with continuous rates."""
    tau = t_expiry - t_start
    v = european_value(1.0, tau, 1.0, interest_rate, dividend_yield, volatility, opt_type.value)
    return stock_price * np.exp(-dividend_yield * t_start) * v


def test_forward_start_matches_rubinstein():
    """A single period reproduces the closed form for a call and a put."""
    t_start = (start_dt - value_dt) / 365.0
    t_expiry = (expiry_dt - value_dt) / 365.0
    for opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
        option = EquityForwardStartOption(start_dt, expiry_dt, opt_type)
        v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
        assert v == pytest.approx(_rubinstein(opt_type, t_start, t_expiry), rel=1e-10)


def test_forward_start_put_call_parity():

    call = EquityForwardStartOption(
        start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL
    )
    put = EquityForwardStartOption(
        start_dt, expiry_dt, OptionTypes.EUROPEAN_PUT
    )

    call_value = call.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    put_value = put.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    dq_start = dividend_curve.df(start_dt)
    dq_end = dividend_curve.df(expiry_dt)

    df_start = discount_curve.df(start_dt)
    df_end = discount_curve.df(expiry_dt)

    expected = stock_price * (
        dq_end - dq_start * df_end / df_start
    )

    assert call_value - put_value == pytest.approx(
        expected, rel=1e-10, abs=1e-10
    )


def test_forward_start_rejects_valuation_after_start():

    option = EquityForwardStartOption(
        start_dt, expiry_dt, OptionTypes.EUROPEAN_CALL
    )

    later_dt = start_dt.add_days(1)

    with pytest.raises(FinError, match="after strike-setting date"):
        option.value(
            later_dt,
            stock_price,
            FlatDiscountCurve(later_dt, interest_rate),
            FlatDiscountCurve(later_dt, dividend_yield),
            model,
        )

    with pytest.raises(FinError, match="Expiry date must be after start date"):
        EquityForwardStartOption(
            expiry_dt,
            start_dt,
            OptionTypes.EUROPEAN_CALL,
        )


def test_forward_start_repr():

    option = EquityForwardStartOption(
        start_dt,
        expiry_dt,
        OptionTypes.EUROPEAN_PUT,
    )

    text = repr(option)

    assert "EquityForwardStartOption" in text
    assert "START_DATE" in text
    assert "FINAL EXPIRY DATE" in text
    assert "OPTION_TYPE" in text


def test_forward_start_on_start_date():

    option = EquityForwardStartOption(
        start_dt,
        expiry_dt,
        OptionTypes.EUROPEAN_CALL,
    )

    dc = FlatDiscountCurve(start_dt, interest_rate)
    qc = FlatDiscountCurve(start_dt, dividend_yield)

    value = option.value(
        start_dt,
        stock_price,
        dc,
        qc,
        model,
    )

    tau = (expiry_dt - start_dt) / 365.0

    expected = european_value(
        stock_price,
        tau,
        stock_price,
        interest_rate,
        dividend_yield,
        volatility,
        OptionTypes.EUROPEAN_CALL.value,
    )

    assert value == pytest.approx(expected, rel=1e-10)
