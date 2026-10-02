"""Floating-lookback limits checked against deterministic path cashflows."""

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.global_types import OptionTypes


@pytest.mark.parametrize("opt_type", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT])
@pytest.mark.parametrize("r,q", [(0.05, 0.01), (0.01, 0.05), (0.03, 0.03)])
def test_zero_volatility_matches_discounted_path_extreme_payoff(opt_type, r, q):
    """Zero volatility yields a known monotone or constant forward path."""
    today = Date(2, 1, 2026)
    expiry = Date(2, 1, 2027)
    times = np.linspace(0.0, 1.0, 1001)
    path = 100.0 * np.exp((r - q) * times)
    historical = 99.0 if opt_type == OptionTypes.EUROPEAN_CALL else 101.0
    if opt_type == OptionTypes.EUROPEAN_CALL:
        payoff = path[-1] - min(historical, np.min(path))
    else:
        payoff = max(historical, np.max(path)) - path[-1]
    option = EquityFloatLookbackOption(expiry, opt_type)
    value = option.value(today, 100.0, FlatDiscountCurve(today, r), FlatDiscountCurve(today, q), BlackScholes(0.0), historical)
    assert value == pytest.approx(np.exp(-r) * payoff, abs=1e-10)


@pytest.mark.parametrize("opt_type,historical", [(OptionTypes.EUROPEAN_CALL, 90.0), (OptionTypes.EUROPEAN_PUT, 110.0)])
def test_expiry_uses_known_historical_extreme(opt_type, historical):
    """The final spot and retained historical extreme determine the payoff."""
    expiry = Date(2, 1, 2027)
    curve = FlatDiscountCurve(expiry, .03)
    option = EquityFloatLookbackOption(expiry, opt_type)
    assert option.value(expiry, 100.0, curve, curve, BlackScholes(.2), historical) == 10.0


def test_post_expiry_fails_before_pricing():
    """Expired contracts cannot be repriced with a negative time interval."""
    expiry = Date(2, 1, 2027)
    after = expiry.add_days(1)
    curve = FlatDiscountCurve(after, .03)
    option = EquityFloatLookbackOption(expiry, OptionTypes.EUROPEAN_CALL)
    with pytest.raises(FinError):
        option.value(after, 100.0, curve, curve, BlackScholes(.2), 90.0)
