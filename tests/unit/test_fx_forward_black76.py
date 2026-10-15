"""Regression tests for FX forward Black-76 valuation (issue #200)."""

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black import Black
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize(
    "option_type", (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT)
)
def test_fx_black76_inferred_forward_matches_direct_black(option_type):
    value_dt = Date(2, 1, 2026)
    expiry_dt = Date(2, 1, 2027)
    spot = 1.12
    strike = 1.10
    domestic_curve = FlatDiscountCurve(value_dt, 0.045)
    foreign_curve = FlatDiscountCurve(value_dt, 0.025)
    option = FXVanillaOption(
        expiry_dt, strike, "EURUSD", option_type, 1_000_000, "USD"
    )
    model = Black(0.18)

    df_dom = domestic_curve.df(expiry_dt)
    df_for = foreign_curve.df(expiry_dt)
    forward = spot * df_for / df_dom
    t_exp = (expiry_dt - value_dt) / 365.0
    expected = model.value(forward, strike, t_exp, df_dom, option_type)

    result = option.value(value_dt, spot, domestic_curve, foreign_curve, model)
    assert result["v"] == pytest.approx(expected, rel=1e-12)
    assert result["pips_dom"] == pytest.approx(expected, rel=1e-12)


def test_fx_black76_accepts_market_forward_quote():
    value_dt, expiry_dt = Date(2, 1, 2026), Date(2, 1, 2027)
    domestic_curve = FlatDiscountCurve(value_dt, 0.04)
    foreign_curve = FlatDiscountCurve(value_dt, 0.02)
    option = FXVanillaOption(
        expiry_dt, 1.15, "EURUSD", OptionTypes.EUROPEAN_CALL, 100_000, "USD"
    )
    model = Black(0.2)
    market_forward = 1.18
    expected = model.value(
        market_forward, 1.15, 1.0,
        domestic_curve.df(expiry_dt), OptionTypes.EUROPEAN_CALL,
    )

    actual = option.value(
        value_dt, 1.13, domestic_curve, foreign_curve, model,
        forward_fx_rate=market_forward,
    )
    assert actual["v"] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("bad_forward", (0.0, -1.0))
def test_fx_black76_rejects_nonpositive_forward(bad_forward):
    value_dt, expiry_dt = Date(2, 1, 2026), Date(2, 1, 2027)
    curve = FlatDiscountCurve(value_dt, 0.03)
    option = FXVanillaOption(
        expiry_dt, 1.0, "EURUSD", OptionTypes.EUROPEAN_CALL, 100, "USD"
    )
    with pytest.raises(FinError, match="forward"):
        option.value(value_dt, 1.0, curve, curve, Black(0.2),
                     forward_fx_rate=bad_forward)


def test_fx_forward_not_silently_ignored_by_spot_model():
    value_dt, expiry_dt = Date(2, 1, 2026), Date(2, 1, 2027)
    curve = FlatDiscountCurve(value_dt, 0.03)
    option = FXVanillaOption(
        expiry_dt, 1.0, "EURUSD", OptionTypes.EUROPEAN_CALL, 100, "USD"
    )
    with pytest.raises(FinError, match="Black model"):
        option.value(value_dt, 1.0, curve, curve, BlackScholes(0.2),
                     forward_fx_rate=1.1)


def test_fx_black76_rejects_american_option():
    value_dt, expiry_dt = Date(2, 1, 2026), Date(2, 1, 2027)
    curve = FlatDiscountCurve(value_dt, 0.03)
    option = FXVanillaOption(
        expiry_dt, 1.0, "EURUSD", OptionTypes.AMERICAN_PUT, 100, "USD"
    )
    with pytest.raises(FinError, match="European"):
        option.value(value_dt, 1.0, curve, curve, Black(0.2))


@pytest.mark.parametrize(
    "option_type", (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT)
)
def test_fx_black76_matches_spot_model_under_interest_parity(option_type):
    value_dt, expiry_dt = Date(2, 1, 2026), Date(2, 1, 2027)
    domestic_curve = FlatDiscountCurve(value_dt, 0.045)
    foreign_curve = FlatDiscountCurve(value_dt, 0.025)
    option = FXVanillaOption(
        expiry_dt, 1.10, "EURUSD", option_type, 100_000, "USD"
    )

    forward_model_price = option.value(
        value_dt, 1.12, domestic_curve, foreign_curve, Black(0.18)
    )["v"]
    spot_model_price = option.value(
        value_dt, 1.12, domestic_curve, foreign_curve, BlackScholes(0.18)
    )["v"]
    assert forward_model_price == pytest.approx(spot_model_price, rel=1e-12)
