import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize("currency", ["USD", "EUR"])
@pytest.mark.parametrize("option_type", [OptionTypes.DIGITAL_CALL, OptionTypes.DIGITAL_PUT])
def test_single_zero_volatility_matches_discounted_forward_event(currency, option_type):
    value_dt = Date(2, 1, 2026)
    expiry_dt = Date(2, 1, 2027)
    domestic = FlatDiscountCurve(value_dt, 0.05)
    foreign = FlatDiscountCurve(value_dt, 0.02)
    spots = np.array([0.9, 1.2, 1.5])
    forward = spots * foreign.df(expiry_dt) / domestic.df(expiry_dt)
    event = forward > 1.2 if option_type == OptionTypes.DIGITAL_CALL else forward < 1.2
    amount = domestic.df(expiry_dt) if currency == "USD" else spots * foreign.df(expiry_dt)
    option = FXDigitalOption(expiry_dt, 1.2, "EURUSD", option_type, 100.0, currency)
    with np.errstate(divide="raise", invalid="raise"):
        actual = option.value(value_dt, spots, domestic, foreign, BlackScholes(0.0))
    np.testing.assert_allclose(actual, 100.0 * amount * event, rtol=1e-12)


@pytest.mark.parametrize("currency", ["USD", "EUR"])
@pytest.mark.parametrize("at_expiry", [False, True])
def test_double_zero_variance_vector_limits(currency, at_expiry):
    value_dt = Date(2, 1, 2026)
    expiry_dt = value_dt if at_expiry else Date(2, 1, 2027)
    domestic = FlatDiscountCurve(value_dt, 0.0)
    foreign = FlatDiscountCurve(value_dt, 0.0)
    spots = np.array([0.9, 1.0, 1.2, 1.4, 1.5])
    option = FXDoubleDigitalOption(expiry_dt, 1.4, 1.0, "EURUSD", 100.0, currency)
    expected = 100.0 * np.array([0.0, 0.5, 1.0, 0.5, 0.0])
    if currency == "EUR":
        expected *= spots
    with np.errstate(divide="raise", invalid="raise"):
        actual = option.value(value_dt, spots, domestic, foreign, BlackScholes(0.0))
    # The library's Hull CDF approximation is documented to six decimals.
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("option_type", [OptionTypes.DIGITAL_CALL, OptionTypes.DIGITAL_PUT])
def test_single_expiry_uses_spot_despite_nonzero_curve_rates(option_type):
    expiry_dt = Date(2, 1, 2026)
    domestic = FlatDiscountCurve(expiry_dt, 0.05)
    foreign = FlatDiscountCurve(expiry_dt, 0.02)
    option = FXDigitalOption(expiry_dt, 1.2, "EURUSD", option_type, 100.0, "USD")
    with np.errstate(divide="raise", invalid="raise"):
        actual = option.value(expiry_dt, 1.2, domestic, foreign, BlackScholes(0.2))
    assert actual == pytest.approx(50.0, abs=1e-6)
