###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinSchedule import FinSchedule
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinDateAdjust():

    startDate = FinDate(28, 2, 2008)
    endDate = FinDate(28, 2, 2011)

    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.NONE
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    testCases.header("NO ADJUSTMENTS", "DATE")
    schedule = FinSchedule(startDate,
                           endDate,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND FOLLOWING", "DATE")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(startDate,
                           endDate,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(startDate,
                           endDate,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

    testCases.banner("")
    testCases.header("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING",
                     "DATE")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.UNITED_STATES
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    startDate = FinDate(4, 7, 2008)
    endDate = FinDate(4, 7, 2011)

    schedule = FinSchedule(startDate,
                           endDate,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print("Date:", dt)

###############################################################################


test_FinDateAdjust()
testCases.compareTestCases()
