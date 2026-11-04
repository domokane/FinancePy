# Copyright (C) 2026 Dominic O'Kane

import numpy as np
import pytest

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.products.rates.swaps.fixed_fixed_cross_ccy import FixedFixedCrossCurrencySwap
from financepy.products.rates.swaps.fixed_fixed_x_ccy_swap import FixedFixedXCcySwap
from financepy.products.rates.swaps.fixed_float_cross_ccy_swap import FixedFloatCrossCurrencySwap
from financepy.products.rates.swaps.ibor_ibor_x_ccy_swap import IborIborXCcySwap
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


SWAP_TYPES = (
    FixedFloatCrossCurrencySwap,
    FixedFixedXCcySwap,
    FixedFixedCrossCurrencySwap,
    IborIborXCcySwap,
)


def _build_curve(value_dt):
    dates = [
        value_dt,
        Date(1, 7, 2024),
        Date(1, 1, 2025),
        Date(1, 7, 2025),
    ]
    dfs = np.array([1.0, 0.99, 0.97, 0.95])
    return DiscountCurve(value_dt, dates, dfs, InterpTypes.FLAT_FWD_RATES)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_cross_currency_swap_value_initializes_discount_factor_before_fixed_leg(
    swap_type,
):
    value_dt = Date(1, 1, 2024)
    curve = _build_curve(value_dt)
    swap = swap_type(
        effective_dt=value_dt,
        term_dt_or_tenor="18M",
        fixed_leg_type=SwapTypes.RECEIVE,
        fixed_cpn=0.02,
        fixed_freq_type=FrequencyTypes.SEMI_ANNUAL,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.ACT_360,
        cal_type=CalendarTypes.WEEKEND,
    )

    value = swap.value(value_dt, curve, curve)

    assert np.isfinite(value)
    assert np.isfinite(swap.fixed_leg_value(value_dt, curve))


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_cross_currency_swap_value_is_stable_and_changes_sign_with_leg_direction(
    swap_type,
):
    value_dt = Date(1, 1, 2024)
    curve = _build_curve(value_dt)
    common_args = dict(
        effective_dt=value_dt,
        term_dt_or_tenor="18M",
        fixed_cpn=0.02,
        fixed_freq_type=FrequencyTypes.SEMI_ANNUAL,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.ACT_360,
        cal_type=CalendarTypes.WEEKEND,
    )

    receive_swap = swap_type(fixed_leg_type=SwapTypes.RECEIVE, **common_args)
    pay_swap = swap_type(fixed_leg_type=SwapTypes.PAY, **common_args)

    receive_value = receive_swap.value(value_dt, curve, curve)
    repeated_value = receive_swap.value(value_dt, curve, curve)
    pay_value = pay_swap.value(value_dt, curve, curve)

    assert repeated_value == pytest.approx(receive_value)
    assert pay_value == pytest.approx(-receive_value)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_cross_currency_swap_matches_ibor_swap_on_a_single_curve(swap_type):
    value_dt = Date(1, 1, 2024)
    curve = _build_curve(value_dt)
    common_args = dict(
        effective_dt=value_dt,
        term_dt_or_tenor="18M",
        fixed_leg_type=SwapTypes.PAY,
        fixed_cpn=0.02,
        fixed_freq_type=FrequencyTypes.SEMI_ANNUAL,
        fixed_dc_type=DayCountTypes.ACT_360,
        float_freq_type=FrequencyTypes.QUARTERLY,
        float_dc_type=DayCountTypes.ACT_360,
        cal_type=CalendarTypes.WEEKEND,
    )
    swap = swap_type(**common_args)
    reference = IborSwap(**common_args)

    assert swap.value(value_dt, curve, curve) == pytest.approx(reference.value(value_dt, curve, curve), rel=1e-10)
    assert swap.swap_rate(value_dt, curve) == pytest.approx(reference.swap_rate(value_dt, curve), rel=1e-10)
