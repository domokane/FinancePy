########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_rainbow_option import (
    EquityRainbowOption,
    EquityRainbowOptionTypes,
)
from financepy.products.fx.fx_double_one_touch_option import FXDoubleOneTouchOption
from financepy.products.fx.fx_rainbow_option import FXRainbowOption, FXRainbowOptionTypes
from financepy.utils.date import Date
from financepy.utils.global_types import DoubleBarrierTypes

value_dt = Date(1, 1, 2026)
expiry_dt = Date(1, 1, 2027)
r_dom = 0.03
domestic_curve = FlatDiscountCurve(value_dt, r_dom)
spot_fx_rates = np.array([100.0, 95.0])
foreign_rates = np.array([0.02, 0.04])
volatilities = np.array([0.30, 0.20])
betas = np.array([0.9, 0.5])
strike = 100.0

########################################################################################


@pytest.mark.parametrize(
    "fx_type, equity_type",
    [
        (FXRainbowOptionTypes.CALL_ON_MAXIMUM, EquityRainbowOptionTypes.CALL_ON_MAXIMUM),
        (FXRainbowOptionTypes.PUT_ON_MAXIMUM, EquityRainbowOptionTypes.PUT_ON_MAXIMUM),
        (FXRainbowOptionTypes.CALL_ON_MINIMUM, EquityRainbowOptionTypes.CALL_ON_MINIMUM),
        (FXRainbowOptionTypes.PUT_ON_MINIMUM, EquityRainbowOptionTypes.PUT_ON_MINIMUM),
    ],
)
def test_fx_rainbow_matches_equity_rainbow(fx_type, equity_type):
    """With the foreign rates as dividend yields and a correlation equal to the
    product of the betas, the FX rainbow is the equity rainbow (Stulz 1982)."""
    fx_option = FXRainbowOption(expiry_dt, fx_type, np.array([strike]), 2)
    v_fx = fx_option.value(
        value_dt, spot_fx_rates, domestic_curve, foreign_rates, volatilities, betas
    )
    rho = betas[0] * betas[1]
    corr_matrix = np.array([[1.0, rho], [rho, 1.0]])
    dividend_curves = [FlatDiscountCurve(value_dt, q) for q in foreign_rates]
    equity_option = EquityRainbowOption(expiry_dt, equity_type, [strike], 2)
    v_equity = equity_option.value(
        value_dt, spot_fx_rates, domestic_curve, dividend_curves, volatilities, corr_matrix
    )
    assert v_fx == pytest.approx(v_equity, rel=1e-10)


@pytest.mark.parametrize(
    "fx_type",
    [FXRainbowOptionTypes.CALL_ON_MAXIMUM, FXRainbowOptionTypes.PUT_ON_MINIMUM],
)
def test_fx_rainbow_analytic_agrees_with_monte_carlo(fx_type):
    """The analytic value and the one-factor Monte Carlo use the same correlation,
    including when the two betas differ."""
    fx_option = FXRainbowOption(expiry_dt, fx_type, np.array([strike]), 2)
    v = fx_option.value(
        value_dt, spot_fx_rates, domestic_curve, foreign_rates, volatilities, betas
    )
    v_mc = fx_option.value_mc(
        value_dt,
        expiry_dt,
        spot_fx_rates,
        domestic_curve,
        foreign_rates,
        volatilities,
        betas,
        200000,
        42,
    )
    assert v == pytest.approx(v_mc, rel=0.02)


def test_fx_option_rho_is_the_domestic_rate_sensitivity():
    """The base-class rho bumps the domestic curve in parallel; it used to call a
    curve method that does not exist."""
    spot_fx_rate, volatility, r_for = 1.10, 0.12, 0.05
    foreign_curve = FlatDiscountCurve(value_dt, r_for)
    model = BlackScholes(volatility)
    option = FXDoubleOneTouchOption(expiry_dt, DoubleBarrierTypes.KNOCK_OUT, 1.00, 1.22, 1.0)
    rho = option.rho(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)
    bump = 1e-4
    v = option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)
    v_up = option.value(
        value_dt, spot_fx_rate, domestic_curve.bump_parallel(bump), foreign_curve, model
    )
    assert rho == pytest.approx((v_up - v) / bump, rel=1e-6)
    assert np.isfinite(rho)
