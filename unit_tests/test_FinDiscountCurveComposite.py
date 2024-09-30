from helpers import buildIborSingleCurve
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_vars import g_basis_point, g_percent
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve_pwf_onf import DiscountCurvePWFONF
from financepy.market.curves.composite_discount_curve import (
    CompositeDiscountCurve,
)
from financepy.products.rates.ibor_swap import IborSwap


def test_composite_discount_curve_can_value_trades():
    """Test that we can use a composite discount curve to value trades"""
    valuation_date = Date(6, 10, 2022)
    base_curve = DiscountCurveFlat(valuation_date, 0.02)

    bump_start_date = Date(6, 10, 2023)
    bump_end_date = Date(6, 10, 2024)
    bump_size = 1.0 * g_percent
    fwd_rate_shock = DiscountCurvePWFONF.brick_wall_curve(
        base_curve.value_dt, bump_start_date, bump_end_date, bump_size
    )
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(valuation_date)
    atm = trade.swap_rate(valuation_date, composite_curve)
    trade.set_fixed_rate(atm)
    value = trade.value(valuation_date, composite_curve)

    expected_atm = 0.02219285487480169
    expected_value = 0.0

    assert abs(atm - expected_atm) < 1e-5
    assert abs(value - expected_value) < 1e-5


def test_zero_bump_has_no_effect_on_base_discount_curve():
    """Test that adding a zero bump to a DiscountCurve does not affect valuation"""
    valuation_date = Date(6, 10, 2022)
    base_curve = DiscountCurveFlat(valuation_date, 0.02)

    bump_start_date = Date(6, 10, 2023)
    bump_end_date = Date(6, 10, 2024)
    bump_size = 0.0 * g_percent
    fwd_rate_shock = DiscountCurvePWFONF.brick_wall_curve(
        base_curve.value_dt, bump_start_date, bump_end_date, bump_size
    )
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(valuation_date)
    atm_base = trade.swap_rate(valuation_date, base_curve)
    atm_comp = trade.swap_rate(valuation_date, composite_curve)

    assert abs(atm_base - atm_comp) < 1e-8

    trade.set_fixed_rate(atm_base)
    value_base = trade.value(valuation_date, base_curve)
    value_comp = trade.value(valuation_date, composite_curve)

    assert abs(value_base - value_comp) < 1e-8


def test_zero_bump_has_no_effect_on_base_ibor_single_curve():
    """Test that adding a zero bump to an IborSingleCurve does not affect valuation"""
    valuation_date = Date(6, 10, 2022)
    base_curve = buildIborSingleCurve(valuation_date, last_tenor="10Y")

    bump_start_date = Date(6, 10, 2023)
    bump_end_date = Date(6, 10, 2024)
    bump_size = 0.0 * g_percent
    fwd_rate_shock = DiscountCurvePWFONF.brick_wall_curve(
        base_curve.value_dt, bump_start_date, bump_end_date, bump_size
    )
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(valuation_date)
    atm_base = trade.swap_rate(valuation_date, base_curve)
    atm_comp = trade.swap_rate(valuation_date, composite_curve)

    assert abs(atm_base - atm_comp) < 1e-8

    trade.set_fixed_rate(atm_base)
    value_base = trade.value(valuation_date, base_curve)
    value_comp = trade.value(valuation_date, composite_curve)

    assert abs(value_base - value_comp) < 1e-8


def _create_test_swap(valuation_date):
    spot_days = 2
    cal = CalendarTypes.UNITED_KINGDOM
    settlement_date = valuation_date.add_weekdays(spot_days)
    swapType = SwapTypes.PAY
    fixedDCCType = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freqType = FrequencyTypes.SEMI_ANNUAL

    trade = IborSwap(
        settlement_date.add_tenor("1Y"),
        "5Y",
        swapType,
        2.0 * g_percent,
        fixed_freqType,
        fixedDCCType,
        cal_type=cal,
        notional=10000,
    )

    return trade


if __name__ == "__main__":
    test_zero_bump_has_no_effect_on_base_ibor_single_curve()
