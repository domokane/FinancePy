###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# TODO Set up test cases correctly

import sys
sys.path.append("..")

from financepy.utils.Calendar import FinDateGenRuleTypes
from financepy.utils.Calendar import FinBusDayAdjustTypes
from financepy.utils.DayCount import FinDayCountTypes
from financepy.utils.Calendar import FinCalendarTypes
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.Date import Date, setDateFormatType, FinDateFormatTypes
from financepy.products.bonds.BondAnnuity import bond_annuity

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBondAnnuity():

    settlement_date = Date(20, 6, 2018)

    #   print("==============================================================")
    #   print("SEMI-ANNUAL FREQUENCY")
    #   print("==============================================================")

    maturity_date = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360
    face = 1000000

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("===============================================================")
#    print("QUARTERLY FREQUENCY")
#    print("===============================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("MONTHLY FREQUENCY")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FinFrequencyTypes.MONTHLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FinFrequencyTypes.ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("BACKWARD GEN WITH SHORT END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FinFrequencyTypes.ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN WITH LONG END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = bond_annuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basisType,
        face)

    annuity.calculateFlowDatesPayments(settlement_date)

    testCases.header("Date", "Flow")
    numFlows = len(annuity._flow_dates)
    for i in range(1, numFlows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

##########################################################################


test_FinBondAnnuity()
testCases.compareTestCases()
