###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.date import Date
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_date_adjust():

    start_date = Date(28, 2, 2008)
    end_date = Date(28, 2, 2011)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.NONE
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    testCases.header("NO ADJUSTMENTS", "DATE")
    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjusted_dates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjusted_dates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjusted_dates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING",
                     "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.UNITED_STATES
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    start_date = Date(4, 7, 2008)
    end_date = Date(4, 7, 2011)

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    for dt in schedule._adjusted_dates:
        testCases.print("Date:", dt)

###############################################################################


test_date_adjust()
testCases.compareTestCases()
