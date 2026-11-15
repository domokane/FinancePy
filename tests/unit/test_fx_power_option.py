# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_european_analytic import european_value
from financepy.products.fx.fx_power_option import FXPowerOption
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes


def test_fx_power_option_power_one_matches_garman_kohlhagen():
    """Verify that power=1.0 exactly reproduces standard Garman-Kohlhagen FX options."""
    value_dt = Date(1, 6, 2026)
    expiry_dt = Date(1, 12, 2026)

    spot_fx_rate = 1.250
    strike_fx_rate = 1.220
    domestic_rate = 0.040
    foreign_rate = 0.025
    volatility = 0.120
    notional = 100000.0

    domestic_curve = FlatDiscountCurve(value_dt, domestic_rate)
    foreign_curve = FlatDiscountCurve(value_dt, foreign_rate)
    model = BlackScholes(volatility)

    t_exp = (expiry_dt - value_dt) / 365.0
    fwd = spot_fx_rate * np.exp((domestic_rate - foreign_rate) * t_exp)

    for opt_type in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        power_option = FXPowerOption(
            expiry_dt,
            strike_fx_rate,
            opt_type,
            power=1.0,
            notional=notional,
            currency_pair="EURUSD",
        )

        power_val = power_option.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        # Closed-form Garman-Kohlhagen
        gk_val = (
            notional
            * np.exp(-domestic_rate * t_exp)
            * european_value(
                fwd, t_exp, strike_fx_rate, 0.0, 0.0, volatility, opt_type.value
            )
        )

        assert np.isclose(power_val, gk_val, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("power", [0.5, 1.0, 1.5, 2.0])
@pytest.mark.parametrize(
    "opt_type", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_fx_power_option_analytical_delta_matches_finite_difference(power, opt_type):
    """Analytical delta must match finite difference bumping within tight relative tolerance."""
    value_dt = Date(1, 3, 2026)
    expiry_dt = Date(1, 9, 2026)

    spot_fx_rate = 1.30
    strike_fx_rate = (spot_fx_rate**power) * 1.03
    domestic_curve = FlatDiscountCurve(value_dt, 0.035)
    foreign_curve = FlatDiscountCurve(value_dt, 0.015)
    model = BlackScholes(0.14)
    notional = 50000.0

    option = FXPowerOption(
        expiry_dt,
        strike_fx_rate,
        opt_type,
        power=power,
        notional=notional,
    )

    h = spot_fx_rate * 1e-6
    v_up = option.value(
        value_dt, spot_fx_rate + h, domestic_curve, foreign_curve, model
    )
    v_dn = option.value(
        value_dt, spot_fx_rate - h, domestic_curve, foreign_curve, model
    )
    bump_delta = (v_up - v_dn) / (2.0 * h)

    analytic_delta = option.delta(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )

    assert analytic_delta == pytest.approx(bump_delta, rel=1e-4, abs=1e-6)


@pytest.mark.parametrize("power", [1.0, 1.5, 2.0])
@pytest.mark.parametrize(
    "opt_type", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_fx_power_option_monte_carlo_agreement(power, opt_type):
    """Analytical price must agree with Monte Carlo simulation within 3 standard errors."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 7, 2026)

    spot_fx_rate = 1.20
    strike_fx_rate = (spot_fx_rate**power) * 1.02
    domestic_curve = FlatDiscountCurve(value_dt, 0.04)
    foreign_curve = FlatDiscountCurve(value_dt, 0.02)
    model = BlackScholes(0.15)
    notional = 1.0

    option = FXPowerOption(
        expiry_dt,
        strike_fx_rate,
        opt_type,
        power=power,
        notional=notional,
    )

    analytic_val = option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )
    mc_result = option.value_mc(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model,
        num_paths=100000,
        seed=12345,
    )

    diff = abs(mc_result.value - analytic_val)
    assert diff <= 3.0 * mc_result.std_err


@pytest.mark.parametrize("power", [0.5, 1.0, 1.5, 2.0])
def test_fx_power_option_put_call_parity(power):
    """Call - Put must satisfy generalized power option put-call parity."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 1, 2027)

    spot = 1.25
    strike = (spot**power) * 1.05
    rd = 0.045
    rf = 0.020
    v = 0.16
    t = 1.0
    notional = 100.0

    dom_curve = FlatDiscountCurve(value_dt, rd)
    for_curve = FlatDiscountCurve(value_dt, rf)
    model = BlackScholes(v)

    call = FXPowerOption(expiry_dt, strike, OptionTypes.EUROPEAN_CALL, power=power, notional=notional)
    put = FXPowerOption(expiry_dt, strike, OptionTypes.EUROPEAN_PUT, power=power, notional=notional)

    call_price = call.value(value_dt, spot, dom_curve, for_curve, model)
    put_price = put.value(value_dt, spot, dom_curve, for_curve, model)

    # Theoretical forward discounted power asset minus discounted strike
    moment_rate = (power - 1.0) * rd - power * rf + 0.5 * power * (power - 1.0) * v * v
    expected_diff = notional * (
        (spot**power) * np.exp(moment_rate * t) - strike * np.exp(-rd * t)
    )

    assert np.isclose(call_price - put_price, expected_diff, rtol=1e-10, atol=1e-10)


def test_fx_power_option_intrinsic_and_zero_vol_limit():
    """As volatility -> 0, option value converges to discounted intrinsic value."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 7, 2026)

    spot = 1.20
    power = 1.5
    strike = (spot**power) * 0.95  # In the money call
    dom_curve = FlatDiscountCurve(value_dt, 0.03)
    for_curve = FlatDiscountCurve(value_dt, 0.01)

    call = FXPowerOption(
        expiry_dt, strike, OptionTypes.EUROPEAN_CALL, power=power, notional=2000.0
    )
    intrinsic = call.intrinsic(value_dt, spot, dom_curve, for_curve)

    # Very low volatility model
    low_vol_model = BlackScholes(1e-6)
    low_vol_val = call.value(value_dt, spot, dom_curve, for_curve, low_vol_model)

    assert np.isclose(low_vol_val, intrinsic, rtol=1e-5, atol=1e-5)


def test_fx_power_option_vectorization():
    """Verify vectorization over spot FX rates for value and delta."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 7, 2026)

    spots = np.array([1.15, 1.20, 1.25, 1.30])
    strike = 1.20
    dom_curve = FlatDiscountCurve(value_dt, 0.04)
    for_curve = FlatDiscountCurve(value_dt, 0.02)
    model = BlackScholes(0.15)

    option = FXPowerOption(
        expiry_dt, strike, OptionTypes.EUROPEAN_CALL, power=1.5, notional=100.0
    )

    vect_vals = option.value(value_dt, spots, dom_curve, for_curve, model)
    vect_deltas = option.delta(value_dt, spots, dom_curve, for_curve, model)

    assert len(vect_vals) == len(spots)
    assert len(vect_deltas) == len(spots)

    for i, s in enumerate(spots):
        scalar_val = option.value(value_dt, s, dom_curve, for_curve, model)
        scalar_delta = option.delta(value_dt, s, dom_curve, for_curve, model)
        assert np.isclose(vect_vals[i], scalar_val, rtol=1e-12, atol=1e-12)
        assert np.isclose(vect_deltas[i], scalar_delta, rtol=1e-12, atol=1e-12)


def test_fx_power_option_notional_scaling():
    """Option price and delta must scale strictly proportionally with notional."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 7, 2026)
    spot = 1.25
    strike = 1.20
    dom_curve = FlatDiscountCurve(value_dt, 0.03)
    for_curve = FlatDiscountCurve(value_dt, 0.01)
    model = BlackScholes(0.12)

    opt1 = FXPowerOption(expiry_dt, strike, OptionTypes.EUROPEAN_CALL, power=1.2, notional=1.0)
    opt100 = FXPowerOption(expiry_dt, strike, OptionTypes.EUROPEAN_CALL, power=1.2, notional=100.0)

    val1 = opt1.value(value_dt, spot, dom_curve, for_curve, model)
    val100 = opt100.value(value_dt, spot, dom_curve, for_curve, model)
    assert np.isclose(val100, val1 * 100.0, rtol=1e-12)

    delta1 = opt1.delta(value_dt, spot, dom_curve, for_curve, model)
    delta100 = opt100.delta(value_dt, spot, dom_curve, for_curve, model)
    assert np.isclose(delta100, delta1 * 100.0, rtol=1e-12)


def test_fx_power_option_validation():
    """Verify input validation and exception raising."""
    expiry_dt = Date(1, 7, 2026)

    # Invalid power (zero, negative, non-finite)
    with pytest.raises(FinError, match="Power must be strictly positive"):
        FXPowerOption(expiry_dt, 1.20, OptionTypes.EUROPEAN_CALL, power=0.0)

    with pytest.raises(FinError, match="Power must be strictly positive"):
        FXPowerOption(expiry_dt, 1.20, OptionTypes.EUROPEAN_CALL, power=-1.5)

    with pytest.raises(FinError, match="Power must be strictly positive"):
        FXPowerOption(expiry_dt, 1.20, OptionTypes.EUROPEAN_CALL, power=np.nan)

    # Invalid strike
    with pytest.raises(FinError):
        FXPowerOption(expiry_dt, -1.20, OptionTypes.EUROPEAN_CALL, power=1.0)

    # Invalid option type
    with pytest.raises(FinError, match="Invalid option type"):
        FXPowerOption(expiry_dt, 1.20, OptionTypes.AMERICAN_CALL, power=1.0)

    # Invalid notional
    with pytest.raises(FinError, match="Notional must be strictly positive"):
        FXPowerOption(expiry_dt, 1.20, OptionTypes.EUROPEAN_CALL, power=1.0, notional=-10.0)

    # Invalid currency pair
    with pytest.raises(FinError, match="Currency pair must be 6 characters"):
        FXPowerOption(
            expiry_dt,
            1.20,
            OptionTypes.EUROPEAN_CALL,
            power=1.0,
            currency_pair="US",
        )


def test_fx_power_option_repr_and_print():
    """Verify string representation and print methods."""
    expiry_dt = Date(1, 7, 2026)
    option = FXPowerOption(
        expiry_dt,
        1.25,
        OptionTypes.EUROPEAN_CALL,
        power=1.5,
        notional=1000.0,
        currency_pair="EURUSD",
    )
    s = repr(option)
    assert "OBJECT TYPE: FXPowerOption" in s
    assert "POWER: 1.5" in s
    assert "CURRENCY PAIR: EURUSD" in s
    assert "NOTIONAL: 1000.0" in s

    # Verify _print does not raise
    option._print()
