# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import add_fp_to_path

from financepy.products.bonds.bond_annuity import BondAnnuity
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes



########################################################################################


def test_bond_annuity():
    """Test"""

    settle_dt = Date(20, 6, 2018)

    #   print("==============================================================")
    #   print("SEMI-ANNUAL FREQUENCY")
    #   print("==============================================================")

    maturity_dt = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    accrual_dc_type = DayCountTypes.ACT_360
    face = 1000000

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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)

    #    print("===============================================================")
    #    print("QUARTERLY FREQUENCY")
    #    print("===============================================================")

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)

    #    print("==================================================================")
    #    print("MONTHLY FREQUENCY")
    #    print("==================================================================")

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.MONTHLY
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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)

    #    print("==================================================================")
    #    print("FORWARD GEN")
    #    print("==================================================================")

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)

    #    print("==================================================================")
    #    print("BACKWARD GEN WITH SHORT END STUB")
    #    print("==================================================================")

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)

    #    print("==================================================================")
    #    print("FORWARD GEN WITH LONG END STUB")
    #    print("==================================================================")

    maturity_dt = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
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

    print("Date", "Flow")
    num_flows = len(annuity.cpn_dts)
    for i in range(1, num_flows):
        dt = annuity.cpn_dts[i]
        flow = annuity.flow_amounts[i]
        print(dt, flow)


########################################################################################

test_bond_annuity()
