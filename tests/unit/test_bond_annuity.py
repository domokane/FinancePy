# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import pytest

from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.products.bonds.bond_annuity import BondAnnuity
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

########################################################################################


def test_semi_annual__bond_annuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    # Semi-Annual Frequency
    maturity_dt = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_dc_type = DayCountTypes.ACT_360
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    annuity = BondAnnuity(
        maturity_dt,
        coupon,
        freq_type,
        accrual_dc_type,
        cal_type,
        bd_type,
        dg_type,
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 2 * 1 + 1
    assert len(annuity.cpn_dts) == 2 * 1 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 25278.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


########################################################################################


def test_quarterly__bond_annuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    # Quarterly Frequency
    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    accrual_dc_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt,
        coupon,
        freq_type,
        accrual_dc_type,
        cal_type,
        bd_type,
        dg_type,
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 4 + 1
    assert len(annuity.cpn_dts) == 10 * 4 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 12778.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


########################################################################################


def test_monthly__bond_annuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    # Monthly Frequency
    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt,
        coupon,
        freq_type,
        basis_type,
        cal_type,
        bd_type,
        dg_type,
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 12 + 1
    assert len(annuity.cpn_dts) == 10 * 12 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 4028.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


########################################################################################


def test_forward_gen__bond_annuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt,
        coupon,
        freq_type,
        basis_type,
        cal_type,
        bd_type,
        dg_type,
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 1 + 1
    assert len(annuity.cpn_dts) == 10 * 1 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert round(annuity.flow_amounts[0]) == 0.0
    assert round(annuity.flow_amounts[-1]) == 50694.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


########################################################################################


def test_forward_gen_with_long_end_stub__bond_annuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt,
        coupon,
        freq_type,
        basis_type,
        cal_type,
        bd_type,
        dg_type,
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 2 + 1
    assert len(annuity.cpn_dts) == 10 * 2 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert round(annuity.flow_amounts[0]) == 0.0
    assert round(annuity.flow_amounts[-1]) == 25417.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


def test_bond_annuity_rebuilds_flows_when_face_changes_at_same_settlement():
    settle_dt = Date(20, 6, 2018)
    annuity = BondAnnuity(
        Date(20, 6, 2019),
        0.05,
        FrequencyTypes.SEMI_ANNUAL,
        DayCountTypes.ACT_360,
    )

    annuity.calculate_payments(settle_dt, 1.0)
    one_unit_coupon = annuity.flow_amounts[1]
    annuity.calculate_payments(settle_dt, 100.0)

    assert annuity.flow_amounts[1] == pytest.approx(100.0 * one_unit_coupon)

    curve = FlatDiscountCurve(settle_dt, 0.0)
    price = annuity.dirty_price_from_discount_curve(settle_dt, curve)
    expected = 100.0 * 0.05 * (183.0 + 182.0) / 360.0
    assert price == pytest.approx(expected)


########################################################################################

test_semi_annual__bond_annuity()
