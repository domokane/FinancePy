from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_vars import G_PERCENT
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.pwf_onf_discount_curve import PWFONFDiscountCurve
from financepy.market.curves.composite_discount_curve import (
    CompositeDiscountCurve,
)
from financepy.products.rates.ibor_swap import IborSwap

from .helpers import build_ibor_single_curve

########################################################################################


def test_composite_discount_curve_can_value_trades():
    """Test that we can use a composite discount curve to value trades"""
    value_dt = Date(6, 10, 2022)
    base_curve = FlatDiscountCurve(value_dt, 0.02)

    bump_start_dt = Date(6, 10, 2023)
    bump_end_dt = Date(6, 10, 2024)
    bump_size = 1.0 * G_PERCENT
    fwd_rate_shock = PWFONFDiscountCurve.brick_wall_curve(base_curve.anchor_dt, bump_start_dt, bump_end_dt, bump_size)
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(value_dt)
    atm = trade.swap_rate(value_dt, composite_curve)
    trade.set_fixed_rate(atm)
    value = trade.value(value_dt, composite_curve)

    expected_atm = 0.02219285487480169
    expected_value = 0.0

    assert abs(atm - expected_atm) < 1e-4
    assert abs(value - expected_value) < 1e-4


########################################################################################


def test_zero_bump_has_no_effect_on_base_discount_curve():
    """Test that adding a zero bump to a DiscountCurve does not affect valuation"""
    value_dt = Date(6, 10, 2022)
    base_curve = FlatDiscountCurve(value_dt, 0.02)

    bump_start_dt = Date(6, 10, 2023)
    bump_end_dt = Date(6, 10, 2024)
    bump_size = 0.0 * G_PERCENT
    fwd_rate_shock = PWFONFDiscountCurve.brick_wall_curve(base_curve.anchor_dt, bump_start_dt, bump_end_dt, bump_size)
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(value_dt)
    atm_base = trade.swap_rate(value_dt, base_curve)
    atm_comp = trade.swap_rate(value_dt, composite_curve)

    assert abs(atm_base - atm_comp) < 1e-8

    trade.set_fixed_rate(atm_base)
    value_base = trade.value(value_dt, base_curve)
    value_comp = trade.value(value_dt, composite_curve)

    assert abs(value_base - value_comp) < 1e-8


########################################################################################


def test_zero_bump_has_no_effect_on_base_ibor_single_curve():
    """Test that adding a zero bump to an IborSingleCurve does not affect valuation"""
    value_dt = Date(6, 10, 2022)
    base_curve = build_ibor_single_curve(value_dt, last_tenor="10Y")

    bump_start_dt = Date(6, 10, 2023)
    bump_end_dt = Date(6, 10, 2024)
    bump_size = 0.0 * G_PERCENT
    fwd_rate_shock = PWFONFDiscountCurve.brick_wall_curve(base_curve.anchor_dt, bump_start_dt, bump_end_dt, bump_size)
    composite_curve = CompositeDiscountCurve([base_curve, fwd_rate_shock])

    trade = _create_test_swap(value_dt)
    atm_base = trade.swap_rate(value_dt, base_curve)
    atm_comp = trade.swap_rate(value_dt, composite_curve)

    assert abs(atm_base - atm_comp) < 1e-8

    trade.set_fixed_rate(atm_base)
    value_base = trade.value(value_dt, base_curve)
    value_comp = trade.value(value_dt, composite_curve)

    assert abs(value_base - value_comp) < 1e-8


########################################################################################


def _create_test_swap(value_dt):

    spot_days = 2
    cal = CalendarTypes.LONDON
    settle_dt = value_dt.add_weekdays(spot_days)
    swap_type = SwapTypes.PAY
    fixed_dcc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq_type = FrequencyTypes.SEMI_ANNUAL

    trade = IborSwap(
        settle_dt.add_tenor("1Y"),
        "5Y",
        swap_type,
        2.0 * G_PERCENT,
        fixed_freq_type,
        fixed_dcc_type,
        cal_type=cal,
        notional=10000,
    )

    return trade


########################################################################################

if __name__ == "__main__":
    test_zero_bump_has_no_effect_on_base_ibor_single_curve()


########################################################################################


def test_composite_df_t_returns_float_for_scalar_and_array_for_array():
    """df_t must have the same return types as the other discount curves"""
    import numpy as np

    value_dt = Date(6, 10, 2022)
    base_curve = FlatDiscountCurve(value_dt, 0.03)
    shock_curve = FlatDiscountCurve(value_dt, 0.01)
    composite_curve = CompositeDiscountCurve([base_curve, shock_curve])

    df_scalar = composite_curve.df_t(1.0)
    assert isinstance(df_scalar, float)
    assert abs(df_scalar - base_curve.df_t(1.0) * shock_curve.df_t(1.0)) < 1e-14

    df_scalar_int = composite_curve.df_t(2)
    assert isinstance(df_scalar_int, float)

    times = np.array([0.5, 1.0, 2.0])
    df_array = composite_curve.df_t(times)
    assert isinstance(df_array, np.ndarray)
    assert df_array.shape == times.shape
    assert np.allclose(df_array, base_curve.df_t(times) * shock_curve.df_t(times))

    df_list = composite_curve.df_t([0.5, 1.0])
    assert isinstance(df_list, np.ndarray)
    assert df_list.shape == (2,)

    # the date interface is unchanged
    df_dt = composite_curve.df(value_dt.add_years(1))
    assert isinstance(df_dt, float)
    assert abs(df_dt - df_scalar) < 1e-12
