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

    settle_date = Date(20, 6, 2018)

    #   print("==============================================================")
    #   print("SEMI-ANNUAL FREQUENCY")
    #   print("==============================================================")

    maturity_date = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360
    face = 1000000

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("===============================================================")
#    print("QUARTERLY FREQUENCY")
#    print("===============================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("MONTHLY FREQUENCY")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("BACKWARD GEN WITH SHORT END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

#    print("==================================================================")
#    print("FORWARD GEN WITH LONG END STUB")
#    print("==================================================================")

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          cal_type,
                          bd_type,
                          dg_type,
                          basis_type)

    annuity.calculate_payments(settle_date, face)

    testCases.header("Date", "Flow")
    num_flows = len(annuity._cpn_dates)
    for i in range(1, num_flows):
        dt = annuity._cpn_dates[i]
        flow = annuity._flow_amounts[i]
        testCases.print(dt, flow)

##########################################################################


test_BondAnnuity()
testCases.compareTestCases()
