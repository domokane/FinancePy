"""Zero-volatility limits using actual calibrated CDS and discount curves."""

import math

import pytest

from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.utils.date import Date


@pytest.mark.parametrize("strike", [0.005, 0.01, 0.02])
def test_adjusted_black_zero_volatility_matches_small_positive_limit(strike):
    """Cover payer/receiver limits on both sides of the adjusted forward spread."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 1, 2027)
    maturity_dt = Date(20, 12, 2030)
    recovery = 0.4
    discount_curve = FlatDiscountCurve(value_dt, 0.03)
    index_curve = CDSCurve(
        value_dt,
        [CDS(value_dt, maturity_dt, 0.01)],
        discount_curve,
        recovery,
    )
    option = CDSIndexOption(expiry_dt, maturity_dt, 0.01, strike)

    zero_vol = option.value_adjusted_black(
        value_dt, index_curve, recovery, discount_curve, 0.0
    )
    small_vol = option.value_adjusted_black(
        value_dt, index_curve, recovery, discount_curve, 1e-7
    )

    assert all(math.isfinite(value) and value >= 0.0 for value in zero_vol)
    assert min(zero_vol) == 0.0
    assert max(zero_vol) > 0.0
    assert zero_vol == pytest.approx(small_vol, rel=1e-10, abs=1e-8)
