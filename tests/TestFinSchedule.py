# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinCalendar import FinDayAdjustTypes, FinDateGenRuleTypes
from financepy.finutils.FinSchedule import FinSchedule
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinSchedule():

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2018, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("SEMI-ANNUAL FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("QUARTERLY FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("MONTHLY FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("FORWARD GEN")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("BACKWARD GEN")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("BACKWARD GEN WITH SHORT END STUB")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("FORWARD GEN WITH LONG END STUB")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    testCases.header("BACKWARD GEN WITH TARGET CALENDAR")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.TARGET
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))


test_FinSchedule()
testCases.compareTestCases()
