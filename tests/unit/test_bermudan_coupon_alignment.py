"""Actual fixed-leg coupons must keep their own payment dates in the tree."""

import numpy as np
import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.bdt_tree import BDTTree
from financepy.models.bk_tree import BKTree
from financepy.models.hw_tree import HWTree
from financepy.products.rates.ibor_bermudan_swaption import IborBermudanSwaption
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import ExerciseTypes, SwapTypes


@pytest.mark.parametrize("tree", [BDTTree(0.01, 50), BKTree(0.01, 0.03, 50), HWTree(0.01, 0.03, 50)])
def test_stub_and_regular_coupons_stay_with_their_payment_dates(tree):
    """The short first period must not receive the final regular coupon."""
    value_dt = Date(2, 1, 2026)
    exercise_dt = Date(17, 3, 2026)
    maturity_dt = Date(30, 9, 2027)
    option = IborBermudanSwaption(
        value_dt, exercise_dt, maturity_dt, SwapTypes.PAY,
        ExerciseTypes.BERMUDAN, 0.04, FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_365F,
    )
    option.value(value_dt, FlatDiscountCurve(value_dt, 0.03), tree)
    leg = option.underlying_swap.fixed_leg
    assert not np.isclose(leg.payments[0], leg.payments[-1])
    np.testing.assert_allclose(option.cpn_flows[1:], np.array(leg.payments) / option.notional)
    expected_times = np.array([(dt - value_dt) / 365.0 for dt in leg.payment_dts])
    np.testing.assert_array_equal(option.cpn_times[1:], expected_times)


def test_european_low_volatility_matches_independent_discounted_cashflows():
    """Near-deterministic receiver value equals fixed less floating cashflows."""
    value_dt = Date(2, 1, 2026)
    exercise_dt = Date(17, 3, 2026)
    maturity_dt = Date(30, 9, 2027)
    curve = FlatDiscountCurve(value_dt, 0.03)
    option = IborBermudanSwaption(
        value_dt, exercise_dt, maturity_dt, SwapTypes.RECEIVE,
        ExerciseTypes.EUROPEAN, 0.08, FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_365F,
    )
    # Daily tree times align exercise and payments exactly, isolating coupon
    # alignment from the tree's date-rounding approximation.
    steps = int(maturity_dt - value_dt)
    value = option.value(value_dt, curve, HWTree(1e-8, 0.03, steps))
    leg = option.underlying_swap.fixed_leg
    fixed_pv = np.dot(leg.payments, curve.df(leg.payment_dts))
    floating_pv = option.notional * (curve.df(exercise_dt) - curve.df(leg.payment_dts[-1]))
    assert value == pytest.approx(max(fixed_pv - floating_pv, 0.0), abs=5.0)
