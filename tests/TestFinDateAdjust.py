###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinSchedule import FinSchedule
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinDateAdjust():

    startDate = FinDate(28, 2, 2008)
    endDate = FinDate(28, 2, 2011)

    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.NONE
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    print("NO ADJUSTMENTS")
    schedule = FinSchedule(startDate,
                           endDate,
                           frequencyType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        print(dt)

    print("")
    print("NO WEEKENDS AND FOLLOWING")
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(startDate,
                           endDate,
                           frequencyType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        print(dt)

    print("")
    print("NO WEEKENDS AND MODIFIED FOLLOWING")
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(startDate,
                           endDate,
                           frequencyType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        print(dt)

    print("")
    print("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING")
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.US
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    startDate = FinDate(4, 7, 2008)
    endDate = FinDate(4, 7, 2011)

    schedule = FinSchedule(startDate,
                           endDate,
                           frequencyType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        print(dt)
###############################################################################


test_FinDateAdjust()
testCases.compareTestCases()
