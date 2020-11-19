###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinSchedule import FinSchedule
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinSchedule():

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2018, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("SEMI-ANNUAL FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("QUARTERLY FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("MONTHLY FREQUENCY")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("FORWARD GEN")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("BACKWARD GEN")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("BACKWARD GEN WITH SHORT END STUB")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    testCases.header("FORWARD GEN WITH LONG END STUB")
    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

    testCases.header("BACKWARD GEN WITH TARGET CALENDAR")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.TARGET
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(
        d1,
        d2,
        freqType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

###############################################################################


test_FinSchedule()
testCases.compareTestCases()
