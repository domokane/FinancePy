###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Date import Date
from financepy.utils.Schedule import Schedule
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.Calendar import FinCalendarTypes
from financepy.utils.Calendar import FinBusDayAdjustTypes
from financepy.utils.Calendar import FinDateGenRuleTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinDateAdjust():

    start_date = Date(28, 2, 2008)
    end_date = Date(28, 2, 2011)

    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.NONE
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    testCases.header("NO ADJUSTMENTS", "DATE")
    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND FOLLOWING", "DATE")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING",
                     "DATE")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.UNITED_STATES
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    start_date = Date(4, 7, 2008)
    end_date = Date(4, 7, 2011)

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

###############################################################################


test_FinDateAdjust()
testCases.compareTestCases()
