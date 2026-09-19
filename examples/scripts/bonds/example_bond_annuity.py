# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.products.bonds.bond_annuity import BondAnnuity
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes

# ============================================================================
# FINANCEPY EXAMPLES - BondAnnuity
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. BOND ANNUITY
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. BOND ANNUITY")
print("=" * 78)

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

