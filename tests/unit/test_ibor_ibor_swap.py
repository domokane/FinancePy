# Copyright (C) 2026 Dominic O'Kane

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.rates.ibor_basis_swap import IborBasisSwap
from financepy.products.rates.swaps.float_float_cross_ccy_swap import FloatFloatCrossCurrencySwap
from financepy.products.rates.swaps.ibor_ibor_swap import IborIborSwap
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


SWAP_TYPES = (IborIborSwap, FloatFloatCrossCurrencySwap)

VALUE_DT = Date(1, 1, 2024)
PAY_FREQ = FrequencyTypes.QUARTERLY
PAY_DC = DayCountTypes.ACT_360
REC_FREQ = FrequencyTypes.SEMI_ANNUAL
REC_DC = DayCountTypes.THIRTY_E_360
SPREAD = 0.001


def _curves():
    discount_curve = FlatDiscountCurve(VALUE_DT, 0.04)
    pay_index_curve = FlatDiscountCurve(VALUE_DT, 0.045)
    rec_index_curve = FlatDiscountCurve(VALUE_DT, 0.05)
    return discount_curve, pay_index_curve, rec_index_curve


def _reference_swap():
    # IborBasisSwap receives leg 1 and pays leg 2, so leg 1 carries the
    # received conventions and the basis spread.
    return IborBasisSwap(
        VALUE_DT,
        "5Y",
        SwapTypes.RECEIVE,
        REC_FREQ,
        REC_DC,
        SPREAD,
        PAY_FREQ,
        PAY_DC,
        0.0,
    )


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_constructs_with_default_calendar(swap_type):
    swap = swap_type(VALUE_DT, "5Y")

    assert swap.maturity_dt == Date(1, 1, 2029)
    assert len(swap._pay_float_dts) == 21
    assert len(swap._rec_float_dts) == 21


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_matches_ibor_basis_swap(swap_type):
    discount_curve, pay_index_curve, rec_index_curve = _curves()
    swap = swap_type(VALUE_DT, "5Y", PAY_FREQ, PAY_DC, REC_FREQ, REC_DC, SPREAD)
    reference = _reference_swap()

    value = swap.value(VALUE_DT, discount_curve, pay_index_curve, rec_index_curve)
    expected = reference.value(VALUE_DT, discount_curve, rec_index_curve, pay_index_curve)

    assert value == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_first_fixings_match_ibor_basis_swap(swap_type):
    discount_curve, pay_index_curve, rec_index_curve = _curves()
    swap = swap_type(VALUE_DT, "5Y", PAY_FREQ, PAY_DC, REC_FREQ, REC_DC, SPREAD)
    reference = _reference_swap()

    value = swap.value(
        VALUE_DT,
        discount_curve,
        pay_index_curve,
        rec_index_curve,
        pay_first_fixing_rate=0.03,
        rec_first_fixing_rate=0.02,
    )
    expected = reference.value(
        VALUE_DT,
        discount_curve,
        rec_index_curve,
        pay_index_curve,
        first_fixing_rate_leg_1=0.02,
        first_fixing_rate_leg_2=0.03,
    )

    assert value == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_identical_legs_value_to_zero(swap_type):
    discount_curve, _, _ = _curves()
    swap = swap_type(VALUE_DT, "5Y")

    value = swap.value(VALUE_DT, discount_curve, discount_curve, discount_curve)

    assert value == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_spread_is_received(swap_type):
    discount_curve, _, _ = _curves()
    swap = swap_type(VALUE_DT, "5Y", basis_swap_spread=SPREAD)

    value = swap.value(VALUE_DT, discount_curve, discount_curve, discount_curve)

    # One basis point of spread on the received leg is worth a positive
    # amount close to the leg PV01.
    assert value > 0.0
    assert value == pytest.approx(SPREAD * swap.notional * 4.5, rel=0.05)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_prints_after_valuation(swap_type, capsys):
    discount_curve, pay_index_curve, rec_index_curve = _curves()
    swap = swap_type(VALUE_DT, "5Y", PAY_FREQ, PAY_DC, REC_FREQ, REC_DC, SPREAD)

    swap.print_float_leg_pv(SwapTypes.PAY)
    assert "Floating Flows not calculated." in capsys.readouterr().out

    swap.value(VALUE_DT, discount_curve, pay_index_curve, rec_index_curve)
    swap.print_float_leg_pv(SwapTypes.PAY)
    swap.print_float_leg_pv(SwapTypes.RECEIVE)
    out = capsys.readouterr().out

    assert "SwapTypes.PAY" in out
    assert "SwapTypes.RECEIVE" in out
    assert "01-JAN-2029" in out
    assert swap_type.__name__ in repr(swap)
