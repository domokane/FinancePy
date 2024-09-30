###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.date import Date


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_dt_adjust():

    start_dt = Date(28, 2, 2008)
    end_dt = Date(28, 2, 2011)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.NONE
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    test_cases.header("NO ADJUSTMENTS", "DATE")
    schedule = Schedule(
        start_dt, end_dt, freq_type, cal_type, bd_type, dg_type
    )

    for dt in schedule.adjusted_dts:
        test_cases.print("Date:", dt)

    test_cases.banner("")
    test_cases.header("NO WEEKENDS AND FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        start_dt, end_dt, freq_type, cal_type, bd_type, dg_type
    )

    for dt in schedule.adjusted_dts:
        test_cases.print("Date:", dt)

    test_cases.banner("")
    test_cases.header("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        start_dt, end_dt, freq_type, cal_type, bd_type, dg_type
    )

    for dt in schedule.adjusted_dts:
        test_cases.print("Date:", dt)

    test_cases.banner("")
    test_cases.header(
        "NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING", "DATE"
    )

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.UNITED_STATES
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    start_dt = Date(4, 7, 2008)
    end_dt = Date(4, 7, 2011)

    schedule = Schedule(
        start_dt, end_dt, freq_type, cal_type, bd_type, dg_type
    )

    for dt in schedule.adjusted_dts:
        test_cases.print("Date:", dt)


###############################################################################


test_dt_adjust()
test_cases.compareTestCases()
