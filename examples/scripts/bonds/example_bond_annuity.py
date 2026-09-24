# ============================================================================
# FINANCEPY EXAMPLES - BondAnnuity
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
#
# This example demonstrates the generation of fixed annuity cash flows using
# different maturities, payment frequencies and date-generation rules.
#
# The examples cover:
#
#   1. Short-dated semi-annual annuity
#   2. Long-dated semi-annual annuity
#   3. Monthly annuity
#   4. Annual annuity with forward date generation
#   5. A repeated annual forward-generation case
#   6. Semi-annual annuity with forward date generation
#
# BondAnnuity generates coupon-style cash flows over a specified schedule.
# The payment amount depends on the coupon rate, face amount, accrual period
# and day-count convention.
# ============================================================================

from financepy.products.bonds.bond_annuity import BondAnnuity
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes

# Common valuation and contract parameters.
settle_dt = Date(20, 6, 2018)
coupon = 0.05
face = 1_000_000

cal_type = CalendarTypes.WEEKEND
bd_type = BusDayAdjustTypes.FOLLOWING
accrual_dc_type = DayCountTypes.ACT_360


# ============================================================================
# 1. SHORT-DATED SEMI-ANNUAL ANNUITY
# ============================================================================
#
# Generate the payment schedule for a one-year annuity with semi-annual
# payments. Dates falling on non-business days are adjusted using the
# Following convention.
# ============================================================================

print("\n" + "=" * 78)
print("1. SHORT-DATED SEMI-ANNUAL ANNUITY")
print("=" * 78)


maturity_dt = Date(20, 6, 2019)
freq_type = FrequencyTypes.SEMI_ANNUAL
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

print("Date", "Flow")

num_flows = len(annuity.cpn_dts)

for i in range(1, num_flows):
    dt = annuity.cpn_dts[i]
    flow = annuity.flow_amounts[i]
    print(dt, flow)


# ============================================================================
# 2. LONG-DATED SEMI-ANNUAL ANNUITY
# ============================================================================
#
# Extend the maturity to ten years while retaining the same semi-annual
# payment frequency and schedule conventions.
#
# This illustrates how BondAnnuity constructs a longer sequence of regular
# coupon-style cash flows.
# ============================================================================

print("\n" + "=" * 78)
print("2. LONG-DATED SEMI-ANNUAL ANNUITY")
print("=" * 78)

maturity_dt = Date(20, 6, 2028)
freq_type = FrequencyTypes.SEMI_ANNUAL
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

print("Date", "Flow")

num_flows = len(annuity.cpn_dts)

for i in range(1, num_flows):
    dt = annuity.cpn_dts[i]
    flow = annuity.flow_amounts[i]
    print(dt, flow)


# ============================================================================
# 3. MONTHLY ANNUITY
# ============================================================================
#
# Generate monthly payments over the same ten-year horizon.
#
# Increasing the payment frequency creates shorter accrual periods and
# therefore smaller individual coupon payments, while producing a much
# larger number of cash flows.
# ============================================================================

print("\n" + "=" * 78)
print("3. MONTHLY ANNUITY")
print("=" * 78)

maturity_dt = Date(20, 6, 2028)
freq_type = FrequencyTypes.MONTHLY
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

print("Date", "Flow")

num_flows = len(annuity.cpn_dts)

for i in range(1, num_flows):
    dt = annuity.cpn_dts[i]
    flow = annuity.flow_amounts[i]
    print(dt, flow)


# ============================================================================
# 4. ANNUAL ANNUITY - FORWARD DATE GENERATION
# ============================================================================
#
# Generate an annual schedule moving forward from the settlement date.
#
# Forward generation is useful when the start of the schedule determines
# the sequence of subsequent payment dates.
# ============================================================================

print("\n" + "=" * 78)
print("4. ANNUAL ANNUITY - FORWARD DATE GENERATION")
print("=" * 78)

maturity_dt = Date(20, 6, 2028)
freq_type = FrequencyTypes.ANNUAL
dg_type = DateGenRuleTypes.FORWARD

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


# ============================================================================
# 5. ANNUAL ANNUITY - REPEATED FORWARD-GENERATION CASE
# ============================================================================
#
# This case is retained from the original example. Its parameters are
# identical to the preceding annual forward-generation example and it
# therefore provides the same schedule and cash flows.
#
# If this section was originally intended to demonstrate a short stub or
# backward date generation, its parameters should be changed accordingly.
# ============================================================================

print("\n" + "=" * 78)
print("5. ANNUAL ANNUITY - REPEATED FORWARD-GENERATION CASE")
print("=" * 78)

maturity_dt = Date(20, 6, 2028)
freq_type = FrequencyTypes.ANNUAL
dg_type = DateGenRuleTypes.FORWARD

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


# ============================================================================
# 6. SEMI-ANNUAL ANNUITY - FORWARD DATE GENERATION
# ============================================================================
#
# Return to semi-annual payments, but generate the schedule forward rather
# than backward. This allows the effect of the date-generation rule to be
# examined while retaining a familiar payment frequency.
# ============================================================================

print("\n" + "=" * 78)
print("6. SEMI-ANNUAL ANNUITY - FORWARD DATE GENERATION")
print("=" * 78)

maturity_dt = Date(20, 6, 2028)
freq_type = FrequencyTypes.SEMI_ANNUAL
dg_type = DateGenRuleTypes.FORWARD

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
