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


def test_fra_pay_fixed_sign_matches_swap_convention():
    """Paying the fixed FRA rate gains when the forward rate is above it, in
    the same way as a pay-fixed IborSwap, and receiving fixed is the mirror."""
    from financepy.products.rates.ibor_swap import IborSwap
    from financepy.utils.frequency import FrequencyTypes
    from financepy.utils.global_types import SwapTypes

    value_dt = Date(2, 1, 2026)
    start_dt = Date(2, 7, 2026)
    maturity_dt = Date(2, 1, 2027)
    curve = FlatDiscountCurve(value_dt, 0.04)
    notional = 1_000_000.0

    alpha = DayCount(DayCountTypes.ACT_360).year_frac(start_dt, maturity_dt)[0]
    fwd = (curve.df(start_dt) / curve.df(maturity_dt) - 1.0) / alpha
    expected = notional * alpha * (fwd - 0.03) * curve.df(maturity_dt) / curve.df(value_dt)
    assert expected > 0.0

    payer = IborFRA(start_dt, maturity_dt, 0.03, DayCountTypes.ACT_360, notional, pay_fixed_rate=True)
    receiver = IborFRA(start_dt, maturity_dt, 0.03, DayCountTypes.ACT_360, notional, pay_fixed_rate=False)

    assert payer.value(value_dt, curve) == pytest.approx(expected, rel=1e-12)
    assert receiver.value(value_dt, curve) == pytest.approx(-expected, rel=1e-12)
    assert payer.valuation_details(value_dt, curve)["value"] == pytest.approx(expected, rel=1e-12)
    assert receiver.valuation_details(value_dt, curve)["value"] == pytest.approx(-expected, rel=1e-12)

    swap = IborSwap(
        value_dt, "5Y", SwapTypes.PAY, 0.03, FrequencyTypes.ANNUAL,
        DayCountTypes.THIRTY_E_360, notional,
    )
    assert swap.value(value_dt, curve) > 0.0
