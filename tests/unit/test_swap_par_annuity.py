import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.products.rates.ois import OIS
from financepy.utils.calendar import CalendarTypes, BusDayAdjustTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


def _swap(cls, leg_type, coupon):
    return cls(Date(2, 1, 2026), Date(2, 1, 2028), leg_type, coupon,
               FrequencyTypes.ANNUAL, DayCountTypes.ACT_365F,
               float_freq_type=FrequencyTypes.ANNUAL,
               float_dc_type=DayCountTypes.ACT_365F,
               cal_type=CalendarTypes.NONE, bd_type=BusDayAdjustTypes.NONE)


@pytest.mark.parametrize("cls", [IborSwap, OIS])
@pytest.mark.parametrize("leg_type", [SwapTypes.PAY, SwapTypes.RECEIVE])
def test_zero_coupon_annuity_and_par_rate_match_independent_cashflows(cls, leg_type):
    value_dt = Date(2, 1, 2026)
    curve = FlatDiscountCurve(value_dt, 0.05)
    swap = _swap(cls, leg_type, 0.0)
    # Two annual ACT/365F coupons on this deliberately unadjusted schedule.
    annuity = curve.df(Date(2, 1, 2027)) + curve.df(Date(2, 1, 2028))
    expected_rate = (1.0 - curve.df(Date(2, 1, 2028))) / annuity
    assert swap.pv01(value_dt, curve) == pytest.approx(annuity)
    assert swap.swap_rate(value_dt, curve) == pytest.approx(expected_rate)
    par_swap = _swap(cls, leg_type, expected_rate)
    assert par_swap.value(value_dt, curve) == pytest.approx(0.0, abs=1e-6)
    if cls is IborSwap:
        assert swap.valuation_details(value_dt, curve)["market_rate"] == pytest.approx(expected_rate)


@pytest.mark.parametrize("leg_type", [SwapTypes.PAY, SwapTypes.RECEIVE])
def test_ois_par_rate_does_not_depend_on_pay_receive_direction(leg_type):
    value_dt = Date(2, 1, 2026)
    curve = FlatDiscountCurve(value_dt, 0.05)
    swap = _swap(OIS, leg_type, 0.04)
    annuity = curve.df(Date(2, 1, 2027)) + curve.df(Date(2, 1, 2028))
    expected = (1.0 - curve.df(Date(2, 1, 2028))) / annuity
    assert swap.swap_rate(value_dt, curve) == pytest.approx(expected)
