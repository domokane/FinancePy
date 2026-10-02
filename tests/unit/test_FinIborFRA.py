"""FRA accrual periods must be positive before any forward-rate division."""

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.rates.ibor_fra import IborFRA
from financepy.utils.date import Date
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.error import FinError


@pytest.mark.parametrize("maturity", [Date(2, 1, 2026), "0D", Date(1, 1, 2026)])
def test_nonpositive_fra_accrual_is_rejected(maturity):
    """Date and tenor constructors reject zero or backwards accrual periods."""
    with pytest.raises(FinError, match="Start date must be before maturity date"):
        IborFRA(Date(2, 1, 2026), maturity, 0.03, DayCountTypes.ACT_360)


def test_positive_fra_at_curve_forward_rate_has_zero_value():
    """A valid FRA still prices at par using the discount-factor ratio."""
    value_dt = Date(2, 1, 2026)
    start_dt = Date(2, 4, 2026)
    maturity_dt = Date(2, 7, 2026)
    curve = FlatDiscountCurve(value_dt, 0.03)
    alpha = DayCount(DayCountTypes.ACT_360).year_frac(start_dt, maturity_dt)[0]
    rate = (curve.df(start_dt) / curve.df(maturity_dt) - 1.0) / alpha
    fra = IborFRA(start_dt, maturity_dt, float(rate), DayCountTypes.ACT_360)
    assert fra.value(value_dt, curve) == pytest.approx(0.0, abs=1e-10)
