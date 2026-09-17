# Copyright (C) 2026 Xamit Kadirbekov

"""Deterministic cash-payment bounds and an independent barrier probability."""

import math

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_double_one_touch_option import FXDoubleOneTouchOption
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import DoubleBarrierTypes


def _prices(spot, lower, upper, days, sigma, domestic_rate, foreign_rate, payment):
    value_dt = Date(1, 1, 2026)
    expiry_dt = value_dt.add_days(days)
    domestic = FlatDiscountCurve(
        value_dt, domestic_rate, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    foreign = FlatDiscountCurve(
        value_dt, foreign_rate, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )
    model = BlackScholes(sigma)
    no_touch = FXDoubleOneTouchOption(
        expiry_dt, DoubleBarrierTypes.KNOCK_OUT, lower, upper, payment
    ).value(value_dt, spot, domestic, foreign, model)
    touch = FXDoubleOneTouchOption(
        expiry_dt, DoubleBarrierTypes.KNOCK_IN, lower, upper, payment
    ).value(value_dt, spot, domestic, foreign, model)
    return no_touch, touch, payment * domestic.df(expiry_dt)


def _assert_cash_bounds_and_parity(no_touch, touch, discounted_payment):
    slack = 1e-12 * discounted_payment
    assert -slack <= no_touch <= discounted_payment + slack
    assert -slack <= touch <= discounted_payment + slack
    assert no_touch + touch == pytest.approx(
        discounted_payment, rel=0, abs=slack
    )


def _normal_interval(lower, upper):
    """Normal probability over an interval without subtracting two near-ones."""
    root_two = math.sqrt(2.0)
    if lower >= 0:
        return 0.5 * (
            math.erfc(lower / root_two) - math.erfc(upper / root_two)
        )
    return 0.5 * (
        math.erfc(-upper / root_two) - math.erfc(-lower / root_two)
    )


def _zero_drift_image_probability(spot, lower, upper, sigma, years):
    """Integrate the absorbing Brownian density using reflected Gaussian images.

    This uses no sine-series coefficients or truncation logic from the pricer.
    For this matrix, 81 images suffice even at its largest diffusion/width ratio.
    Domestic rate = foreign rate + sigma**2/2 makes log-FX drift zero.
    """
    width = math.log(upper / lower)
    position = math.log(spot / lower)
    scale = sigma * math.sqrt(years)
    return math.fsum(
        _normal_interval(
            (-position + 2 * k * width) / scale,
            (width - position + 2 * k * width) / scale,
        )
        - _normal_interval(
            (position + 2 * k * width) / scale,
            (width + position + 2 * k * width) / scale,
        )
        for k in range(-40, 41)
    )


def test_double_one_touch_reported_negative_price():
    """A complementary nonnegative payoff must not have a negative price."""
    sigma = 0.2
    no_touch, touch, discounted_payment = _prices(
        100.0, 80.0, 125.0, 30, sigma, 0.5 * sigma * sigma, 0.0, 1.0
    )
    _assert_cash_bounds_and_parity(no_touch, touch, discounted_payment)
    assert no_touch == pytest.approx(0.9981587590142639, rel=0, abs=2e-10)
    assert touch == pytest.approx(0.00019875572704133946, rel=0, abs=2e-10)


@pytest.mark.parametrize("sigma", [0.1, 0.2, 0.4])
@pytest.mark.parametrize("half_log_width", [0.1, 0.5])
@pytest.mark.parametrize("days", [7, 730])
@pytest.mark.parametrize("position", [0.25, 0.5, 0.75])
@pytest.mark.parametrize("payment", [1.0, 100.0])
def test_double_no_touch_matches_reflected_density(
    sigma, half_log_width, days, position, payment
):
    """Cover zero coefficients, midpoint sine zeros and short-maturity tails."""
    lower = 100.0
    upper = lower * math.exp(2 * half_log_width)
    spot = lower * math.exp(2 * half_log_width * position)
    foreign_rate = 0.03
    domestic_rate = foreign_rate + 0.5 * sigma * sigma
    no_touch, touch, discounted_payment = _prices(
        spot, lower, upper, days, sigma, domestic_rate, foreign_rate, payment
    )
    probability = _zero_drift_image_probability(
        spot, lower, upper, sigma, days / 365.0
    )
    _assert_cash_bounds_and_parity(no_touch, touch, discounted_payment)
    assert no_touch / discounted_payment == pytest.approx(
        probability, rel=0, abs=2e-10
    )


@pytest.mark.parametrize("log_drift", [-0.015, 0.015])
@pytest.mark.parametrize("sigma", [0.1, 0.4])
@pytest.mark.parametrize("days", [30, 730])
@pytest.mark.parametrize("position", [math.sqrt(2) / 4, 0.643])
def test_double_one_touch_nonzero_drift_bounds_and_parity(
    log_drift, sigma, days, position
):
    """Preserve the complementary payoffs away from the zero-drift trigger."""
    no_touch, touch, discounted_payment = _prices(
        100.0 * math.exp(position),
        100.0,
        100.0 * math.exp(1.0),
        days,
        sigma,
        0.5 * sigma * sigma + log_drift,
        0.0,
        100.0,
    )
    _assert_cash_bounds_and_parity(no_touch, touch, discounted_payment)
