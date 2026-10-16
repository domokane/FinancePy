"""Regression coverage for analytical deltas of the existing power option."""
import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_power_option import power_option_delta, power_option_value
from financepy.products.equity.equity_power_option import EquityPowerOption
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize(
    "power,spot,strike",
    [(0.5, 100.0, 10.2), (1.0, 100.0, 105.0),
     (2.0, 10.0, 105.0), (3.0, 5.0, 140.0)],
)
@pytest.mark.parametrize(
    "kind", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_analytic_delta_matches_existing_power_price_finite_difference(
    power, spot, strike, kind
):
    t, rate, dividend, volatility = 1.4, 0.045, 0.012, 0.25
    h = spot * 1e-5

    def price(s):
        return power_option_value(
            s, t, strike, rate, dividend, volatility, kind.value, power
        )

    # The existing pricer uses a ~6-decimal normal-CDF approximation.
    # Price finite differences therefore need a small tolerance.
    bump_delta = (price(spot + h) - price(spot - h)) / (2.0 * h)
    analytic_delta = power_option_delta(
        spot, t, strike, rate, dividend, volatility, kind.value, power
    )
    assert analytic_delta == pytest.approx(bump_delta, rel=5e-5, abs=2e-6)


@pytest.mark.parametrize(
    "kind", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_power_one_delta_matches_black_scholes(kind):
    from scipy.special import ndtr

    s, k, t, r, q, v = 100.0, 105.0, 1.25, 0.04, 0.01, 0.22
    d1 = (np.log(s/k) + (r - q + v*v/2.0)*t) / (v*np.sqrt(t))
    expected = np.exp(-q*t) * (
        ndtr(d1) if kind == OptionTypes.EUROPEAN_CALL else -ndtr(-d1)
    )
    actual = power_option_delta(s, t, k, r, q, v, kind.value, 1.0)
    assert actual == pytest.approx(expected, abs=2e-7)


@pytest.mark.parametrize("power", [0.5, 1.0, 2.0])
@pytest.mark.parametrize(
    "kind", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_equity_power_delta_matches_product_valuation_and_quantity(power, kind):
    today, expiry = Date(1, 1, 2026), Date(1, 7, 2027)
    spot = 10.0
    strike = spot**power * 1.05
    product = EquityPowerOption(expiry, strike, kind, power, num_options=3.0)
    discount = FlatDiscountCurve(today, 0.045)
    dividend = FlatDiscountCurve(today, 0.015)
    model = BlackScholes(0.25)

    h = spot * 1e-5
    up = product.value(today, spot+h, discount, dividend, model)
    dn = product.value(today, spot-h, discount, dividend, model)
    bump_delta = (up-dn)/(2.0*h)

    assert product.delta(today, spot, discount, dividend, model) == pytest.approx(
        bump_delta, rel=5e-5, abs=2e-6
    )


@pytest.mark.parametrize("power", [0.0, -1.0, -0.5, np.nan, np.inf])
def test_equity_power_requires_strictly_positive_finite_power(power):
    with pytest.raises(FinError, match="Power must be strictly positive"):
        EquityPowerOption(
            Date(1, 7, 2027), 100.0, OptionTypes.EUROPEAN_CALL, float(power)
        )


@pytest.mark.parametrize("power", [0.0, -1.0])
def test_model_maintains_nonpositive_power_rejection(power):
    assert np.isnan(power_option_delta(
        100.0, 1.0, 100.0, 0.04, 0.01, 0.2,
        OptionTypes.EUROPEAN_CALL.value, power
    ))


@pytest.mark.parametrize(
    "kind", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]
)
def test_zero_volatility_delta_is_deterministic_away_from_kink(kind):
    spot, t, k, n = 12.0, 1.0, 100.0, 2.0
    analytic = power_option_delta(
        spot, t, k, 0.0, 0.0, 0.0, kind.value, n
    )
    expected = n * spot**(n-1) if kind == OptionTypes.EUROPEAN_CALL else 0.0
    assert analytic == pytest.approx(expected, abs=1e-8)
