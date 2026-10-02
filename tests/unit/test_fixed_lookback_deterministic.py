"""Fixed-lookback limits checked against known deterministic path payoffs."""

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_fixed_lookback_option import EquityFixedLookbackOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize("opt_type", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT])
@pytest.mark.parametrize("r,q", [(0.05, 0.01), (0.01, 0.05), (0.03, 0.03)])
def test_zero_volatility_matches_discounted_path_extreme(opt_type, r, q):
    """The strike payoff uses the historical and deterministic future path."""
    today = Date(2, 1, 2026)
    expiry = Date(2, 1, 2027)
    path = 100.0 * np.exp((r - q) * np.linspace(0.0, 1.0, 1001))
    historical = 101.0 if opt_type == OptionTypes.EUROPEAN_CALL else 99.0
    if opt_type == OptionTypes.EUROPEAN_CALL:
        payoff = max(max(historical, np.max(path)) - 100.0, 0.0)
    else:
        payoff = max(100.0 - min(historical, np.min(path)), 0.0)
    option = EquityFixedLookbackOption(expiry, opt_type, 100.0)
    value = option.value(today, 100.0, FlatDiscountCurve(today, r), FlatDiscountCurve(today, q), BlackScholes(0.0), historical)
    assert value == pytest.approx(np.exp(-r) * payoff, abs=1e-10)


@pytest.mark.parametrize("opt_type,strike,historical,expected", [
    (OptionTypes.EUROPEAN_CALL, 100.0, 110.0, 10.0),
    (OptionTypes.EUROPEAN_CALL, 115.0, 110.0, 0.0),
    (OptionTypes.EUROPEAN_PUT, 100.0, 90.0, 10.0),
    (OptionTypes.EUROPEAN_PUT, 85.0, 90.0, 0.0),
])
def test_expiry_uses_observed_extreme_and_strike(opt_type, strike, historical, expected):
    """Known expiry extrema give the exact in/out-of-money strike payoff."""
    expiry = Date(2, 1, 2027)
    curve = FlatDiscountCurve(expiry, .03)
    option = EquityFixedLookbackOption(expiry, opt_type, strike)
    assert option.value(expiry, 100.0, curve, curve, BlackScholes(.2), historical) == expected
