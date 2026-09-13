# Copyright (C) 2026 Xamit Kadirbekov

from math import erf, exp, log, sqrt

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR


VALUE_DT = Date(1, 1, 2024)
EXPIRY_DT = Date(1, 1, 2025)
SPOT = 1.20
STRIKE = 1.25
LOWER_STRIKE = 1.10
UPPER_STRIKE = 1.40
DOMESTIC_RATE = 0.05
FOREIGN_RATE = 0.01
VOLATILITY = 0.20
TIME = (EXPIRY_DT - VALUE_DT) / G_DAYS_IN_YEAR


def normal_cdf(value):
    return (1.0 + erf(value / sqrt(2.0))) / 2.0


def d1_d2(strike):
    volatility_horizon = VOLATILITY * sqrt(TIME)
    d1 = (
        log(SPOT / strike)
        + (DOMESTIC_RATE - FOREIGN_RATE + VOLATILITY**2 / 2.0) * TIME
    ) / volatility_horizon
    return d1, d1 - volatility_horizon


DOMESTIC_CURVE = FlatDiscountCurve(VALUE_DT, DOMESTIC_RATE)
FOREIGN_CURVE = FlatDiscountCurve(VALUE_DT, FOREIGN_RATE)
MODEL = BlackScholes(VOLATILITY)


@pytest.mark.parametrize(
    "option_type,sign",
    [
        (OptionTypes.DIGITAL_CALL, 1.0),
        (OptionTypes.DIGITAL_PUT, -1.0),
    ],
)
def test_foreign_currency_digital_is_asset_or_nothing(option_type, sign):
    option = FXDigitalOption(
        EXPIRY_DT,
        STRIKE,
        "EURUSD",
        option_type,
        1.0,
        "EUR",
    )

    value = option.value(
        VALUE_DT,
        SPOT,
        DOMESTIC_CURVE,
        FOREIGN_CURVE,
        MODEL,
    )
    d1, _ = d1_d2(STRIKE)
    expected = SPOT * exp(-FOREIGN_RATE * TIME) * normal_cdf(sign * d1)

    # FinancePy's vector normal CDF is a polynomial approximation; the
    # independent oracle uses math.erf.
    assert value == pytest.approx(expected, rel=0.0, abs=2e-7)


@pytest.mark.parametrize(
    "option_type,sign",
    [
        (OptionTypes.DIGITAL_CALL, 1.0),
        (OptionTypes.DIGITAL_PUT, -1.0),
    ],
)
def test_domestic_currency_digital_is_cash_or_nothing(option_type, sign):
    option = FXDigitalOption(
        EXPIRY_DT,
        STRIKE,
        "EURUSD",
        option_type,
        1.0,
        "USD",
    )

    value = option.value(
        VALUE_DT,
        SPOT,
        DOMESTIC_CURVE,
        FOREIGN_CURVE,
        MODEL,
    )
    _, d2 = d1_d2(STRIKE)
    expected = exp(-DOMESTIC_RATE * TIME) * normal_cdf(sign * d2)

    assert value == pytest.approx(expected, rel=0.0, abs=2e-7)


def test_domestic_currency_double_digital_uses_domestic_discounting():
    option = FXDoubleDigitalOption(
        EXPIRY_DT,
        UPPER_STRIKE,
        LOWER_STRIKE,
        "EURUSD",
        1.0,
        "USD",
    )

    value = option.value(
        VALUE_DT,
        SPOT,
        DOMESTIC_CURVE,
        FOREIGN_CURVE,
        MODEL,
    )
    _, lower_d2 = d1_d2(LOWER_STRIKE)
    _, upper_d2 = d1_d2(UPPER_STRIKE)
    expected = exp(-DOMESTIC_RATE * TIME) * (
        normal_cdf(lower_d2) - normal_cdf(upper_d2)
    )

    assert value == pytest.approx(expected, rel=0.0, abs=2e-7)


def test_foreign_currency_double_digital_is_asset_interval():
    option = FXDoubleDigitalOption(
        EXPIRY_DT,
        UPPER_STRIKE,
        LOWER_STRIKE,
        "EURUSD",
        1.0,
        "EUR",
    )

    value = option.value(
        VALUE_DT,
        SPOT,
        DOMESTIC_CURVE,
        FOREIGN_CURVE,
        MODEL,
    )
    lower_d1, _ = d1_d2(LOWER_STRIKE)
    upper_d1, _ = d1_d2(UPPER_STRIKE)
    expected = SPOT * exp(-FOREIGN_RATE * TIME) * (
        normal_cdf(lower_d1) - normal_cdf(upper_d1)
    )

    assert value == pytest.approx(expected, rel=0.0, abs=2e-7)
