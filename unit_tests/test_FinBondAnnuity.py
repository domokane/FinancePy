########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys

sys.path.append("..")

from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.products.bonds.bond_annuity import BondAnnuity


def test_SemiAnnual_BondAnnuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    # Semi-Annual Frequency
    maturity_dt = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt, coupon, freq_type, cal_type, bd_type, dg_type, basis_type
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 2 * 1 + 1
    assert len(annuity.cpn_dts) == 2 * 1 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 25278.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


def test_Quarterly_BondAnnuity():

    settle_dt = Date(20, 6, 2018)
    face = 1000000

    # Quarterly Frequency
    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_dt, coupon, freq_type, cal_type, bd_type, dg_type, basis_type
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 4 + 1
    assert len(annuity.cpn_dts) == 10 * 4 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 12778.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


def test_Monthly_BondAnnuity():

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
        maturity_dt, coupon, freq_type, cal_type, bd_type, dg_type, basis_type
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 12 + 1
    assert len(annuity.cpn_dts) == 10 * 12 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert annuity.flow_amounts[0] == 0.0
    assert round(annuity.flow_amounts[-1]) == 4028.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


def test_ForwardGen_BondAnnuity():

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
        maturity_dt, coupon, freq_type, cal_type, bd_type, dg_type, basis_type
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 1 + 1
    assert len(annuity.cpn_dts) == 10 * 1 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert round(annuity.flow_amounts[0]) == 0.0
    assert round(annuity.flow_amounts[-1]) == 50694.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


def test_ForwardGenWithLongEndStub_BondAnnuity():

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
        maturity_dt, coupon, freq_type, cal_type, bd_type, dg_type, basis_type
    )

    annuity.calculate_payments(settle_dt, face)

    assert len(annuity.flow_amounts) == 10 * 2 + 1
    assert len(annuity.cpn_dts) == 10 * 2 + 1

    assert annuity.cpn_dts[0] == settle_dt
    assert annuity.cpn_dts[-1] == maturity_dt

    assert round(annuity.flow_amounts[0]) == 0.0
    assert round(annuity.flow_amounts[-1]) == 25417.0

    assert annuity.accrued_interest(settle_dt, face) == 0.0


test_SemiAnnual_BondAnnuity()
