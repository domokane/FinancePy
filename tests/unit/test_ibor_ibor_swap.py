# Copyright (C) 2026 Dominic O'Kane

import pytest

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.rates.ibor_basis_swap import IborBasisSwap
from financepy.products.rates.swaps.float_float_cross_ccy_swap import FloatFloatCrossCurrencySwap
from financepy.products.rates.swaps.ibor_ibor_swap import IborIborSwap
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.error import FinError
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes


SWAP_TYPES = (IborIborSwap, FloatFloatCrossCurrencySwap)

EFFECTIVE_DT = Date(1, 1, 2024)
PAY_FREQ = FrequencyTypes.QUARTERLY
PAY_DC = DayCountTypes.ACT_360
REC_FREQ = FrequencyTypes.SEMI_ANNUAL
REC_DC = DayCountTypes.THIRTY_E_360
SPREAD = 0.001
PAY_FIXING = 0.03
REC_FIXING = 0.02

VALUE_DTS = (
    Date(1, 12, 2023),  # forward starting
    Date(1, 1, 2024),  # effective date
    Date(15, 2, 2024),  # seasoned, inside the first period
    Date(1, 7, 2024),  # a payment date on both legs
    Date(1, 1, 2029),  # maturity
    Date(2, 1, 2029),  # after maturity
)


def _curves(value_dt):
    discount_curve = FlatDiscountCurve(value_dt, 0.04)
    pay_index_curve = FlatDiscountCurve(value_dt, 0.045)
    rec_index_curve = FlatDiscountCurve(value_dt, 0.05)
    return discount_curve, pay_index_curve, rec_index_curve


def _swap(swap_type):
    return swap_type(EFFECTIVE_DT, "5Y", PAY_FREQ, PAY_DC, REC_FREQ, REC_DC, SPREAD)


def _reference_swap():
    # IborBasisSwap receives leg 1 and pays leg 2, so leg 1 carries the
    # received conventions and the basis spread.
    return IborBasisSwap(
        EFFECTIVE_DT,
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
def test_ibor_ibor_swap_generates_a_schedule_per_leg(swap_type):
    swap = _swap(swap_type)

    assert swap.maturity_dt == Date(1, 1, 2029)
    assert len(swap._pay_float_dts) == 21
    assert len(swap._rec_float_dts) == 11


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_matches_ibor_basis_swap(swap_type):
    discount_curve, pay_index_curve, rec_index_curve = _curves(EFFECTIVE_DT)

    value = _swap(swap_type).value(EFFECTIVE_DT, discount_curve, pay_index_curve, rec_index_curve)
    expected = _reference_swap().value(EFFECTIVE_DT, discount_curve, rec_index_curve, pay_index_curve)

    assert value == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("value_dt", VALUE_DTS)
@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_with_first_fixings_matches_ibor_basis_swap(swap_type, value_dt):
    discount_curve, pay_index_curve, rec_index_curve = _curves(value_dt)

    value = _swap(swap_type).value(
        value_dt,
        discount_curve,
        pay_index_curve,
        rec_index_curve,
        pay_first_fixing_rate=PAY_FIXING,
        rec_first_fixing_rate=REC_FIXING,
    )
    expected = _reference_swap().value(
        value_dt,
        discount_curve,
        rec_index_curve,
        pay_index_curve,
        first_fixing_rate_leg_1=REC_FIXING,
        first_fixing_rate_leg_2=PAY_FIXING,
    )

    assert value == pytest.approx(expected, rel=1e-10, abs=1e-8)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_values_to_zero_after_maturity(swap_type):
    value_dt = Date(2, 1, 2029)
    discount_curve, pay_index_curve, rec_index_curve = _curves(value_dt)

    assert _swap(swap_type).value(value_dt, discount_curve, pay_index_curve, rec_index_curve) == 0.0


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_identical_legs_value_to_zero(swap_type):
    discount_curve, _, _ = _curves(EFFECTIVE_DT)
    swap = swap_type(EFFECTIVE_DT, "5Y")

    value = swap.value(EFFECTIVE_DT, discount_curve, discount_curve, discount_curve)

    assert value == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_spread_is_received(swap_type):
    discount_curve, _, _ = _curves(EFFECTIVE_DT)
    swap = swap_type(EFFECTIVE_DT, "5Y", basis_swap_spread=SPREAD)

    value = swap.value(EFFECTIVE_DT, discount_curve, discount_curve, discount_curve)

    # Ten basis points of spread on the received leg is worth roughly the
    # spread times the notional times the annuity of a five year swap.
    assert value > 0.0
    assert value == pytest.approx(SPREAD * swap.notional * 4.5, rel=0.05)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_rejects_curves_with_different_anchor_dates(swap_type):
    discount_curve, pay_index_curve, _ = _curves(EFFECTIVE_DT)
    stale_curve = FlatDiscountCurve(Date(1, 12, 2023), 0.05)

    with pytest.raises(FinError):
        _swap(swap_type).value(EFFECTIVE_DT, discount_curve, pay_index_curve, stale_curve)


@pytest.mark.parametrize("swap_type", SWAP_TYPES)
def test_ibor_ibor_swap_prints_each_leg_after_valuation(swap_type, capsys):
    discount_curve, pay_index_curve, rec_index_curve = _curves(EFFECTIVE_DT)
    swap = _swap(swap_type)

    swap.print_float_leg_pv(SwapTypes.PAY)
    assert "Floating Flows not calculated." in capsys.readouterr().out

    swap.value(EFFECTIVE_DT, discount_curve, pay_index_curve, rec_index_curve)

    swap.print_float_leg_pv(SwapTypes.PAY)
    pay_out = capsys.readouterr().out
    swap.print_float_leg_pv(SwapTypes.RECEIVE)
    rec_out = capsys.readouterr().out

    assert "Floating Flows not calculated." not in pay_out
    assert "Floating Flows not calculated." not in rec_out
    assert "FIRST FLOATING RATE PAYMENT IS IMPLIED" in pay_out
    assert "01-APR-2024" in pay_out  # quarterly pay leg
    assert "01-APR-2024" not in rec_out  # semi-annual receive leg
    assert "01-JUL-2024" in rec_out
    assert "01-JAN-2029" in rec_out

    text = repr(swap)
    assert swap_type.__name__ in text
    assert "BASIS SWAP SPREAD" in text
    assert "REC FREQUENCY" in text
