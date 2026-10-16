"""The curve interpolation and floating-coupon accrual clocks are distinct."""
import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.rates.swap_float_leg import SwapFloatLeg
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


@pytest.mark.parametrize("day_count", [
    DayCountTypes.ACT_360, DayCountTypes.THIRTY_E_360,
    DayCountTypes.ACT_365F,
])
def test_float_coupon_forward_uses_index_accrual_clock(day_count):
    effective = Date(5, 1, 2026)
    leg = SwapFloatLeg(
        effective, "6M", SwapTypes.RECEIVE, 0.0,
        FrequencyTypes.QUARTERLY, day_count,
        notional=1_000_000,
    )
    discount = FlatDiscountCurve(effective, .025)
    index = FlatDiscountCurve(effective, .047)
    price = leg.value(effective, discount, index)

    assert len(leg.rates) == 2
    expected = 0.0
    for i, (start, end, pay_dt) in enumerate(zip(
        leg.start_accrued_dts, leg.end_accrued_dts, leg.payment_dts,
    )):
        forward = index.fwd_rate(start, end, accrual_dc_type=day_count)
        assert leg.rates[i] == pytest.approx(forward, rel=1e-12)
        alpha = leg.year_fracs[i]
        expected += 1_000_000 * alpha * forward * discount.df(pay_dt)

    assert price == pytest.approx(expected, rel=1e-11)


def test_float_coupon_first_fixing_still_overrides_only_first_forward():
    value_date = Date(5, 1, 2026)
    leg = SwapFloatLeg(
        value_date, "6M", SwapTypes.PAY, .001,
        FrequencyTypes.QUARTERLY, DayCountTypes.THIRTY_E_360,
        notional=1_000_000,
    )
    curve = FlatDiscountCurve(value_date, .04)
    price = leg.value(value_date, curve, curve, first_fixing_rate=.032)
    assert leg.rates[0] == .032
    assert leg.rates[1] == pytest.approx(
        curve.fwd_rate(leg.start_accrued_dts[1], leg.end_accrued_dts[1],
                       accrual_dc_type=DayCountTypes.THIRTY_E_360)
    )
    assert np.isfinite(price) and price < 0.0

