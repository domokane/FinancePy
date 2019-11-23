# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""
# TODO Set up test cases correctly

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.FinAnnuity import FinAnnuity
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinCalendar import FinDayAdjustTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinAnnuity():

#    print("==================================================================")
#    print("SEMI-ANNUAL FREQUENCY")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2018, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    basisType = FinDayCountTypes.ACT_360

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("QUARTERLY FREQUENCY")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.QUARTERLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("MONTHLY FREQUENCY")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("FOWARD GEN")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("BACKWARD GEN")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("BACKWARD GEN WITH SHORT END STUB")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("FORWARD GEN WITH LONG END STUB")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

#    print("==================================================================")
#    print("BACKWARD GEN WITH TARGET CALENDAR")
#    print("==================================================================")

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2028, 6, 20)
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.TARGET
    busDayAdjustType = FinDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    annuity = FinAnnuity(
        d1,
        d2,
        frequencyType,
        calendarType,
        busDayAdjustType,
        dateGenRuleType,
        basisType)
#    annuity.dump()

##########################################################################


test_FinAnnuity()
testCases.compareTestCases()
