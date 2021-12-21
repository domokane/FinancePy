###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# TODO Set up test cases correctly
import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.bond_annuity import BondAnnuity
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_BondAnnuity():

    settlement_date = Date(20, 6, 2018)

    #   print("==============================================================")
    #   print("SEMI-ANNUAL FREQUENCY")
    #   print("==============================================================")

    maturity_date = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360
    face = 1000000

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type,
                          face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("===============================================================")
#    print("QUARTERLY FREQUENCY")
#    print("===============================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basis_type,
        face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("MONTHLY FREQUENCY")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.MONTHLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type,
                          face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type,
                          face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("BACKWARD GEN WITH SHORT END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type,
                          face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN WITH LONG END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type,
                          face)

    annuity.calculate_payments(settlement_date)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._flow_dates)
    for i in range(1, num_flows):
        dt = annuity._flow_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

##########################################################################


test_BondAnnuity()
testCases.compareTestCases()
