###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinCalendar import FinBusDayAdjustTypes
from financepy.finutils.FinCalendar import FinDateGenRuleTypes
from financepy.finutils.FinSchedule import FinSchedule
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes, FinCalendar
from financepy.finutils.FinDate import FinDate, setDateFormatType, FinDateFormatTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

setDateFormatType(FinDateFormatTypes.UK_LONGEST)

def test_FinSchedule():

    d1 = FinDate(2018, 6, 20)
    d2 = FinDate(2018, 6, 20)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
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

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    for dt in schedule._adjustedDates:
        testCases.print(str(dt))

###############################################################################

def test_FinScheduleAligment(eomFlag):
        
    valuationDate = FinDate(m=3,d=29,y=2005)
    effDate = valuationDate.addTenor("2d")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    calendarType = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = False

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")

#    print(matDate1)
#    print(matDate2)

    myCal = FinCalendar(calendarType)

    adjustedMatDate1 = myCal.adjust(matDate1, busDayAdjustType)
    adjustedMatDate2 = myCal.adjust(matDate2, busDayAdjustType)

#    print(adjustedMatDate1)
#    print(adjustedMatDate2)
    
    sched1 = FinSchedule(effDate,
                         adjustedMatDate1,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)
    
#    print(sched1)
    
    sched2 = FinSchedule(effDate,
                         adjustedMatDate2,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == eomFlag)

###############################################################################

test_FinSchedule()
test_FinScheduleAligment(True)
test_FinScheduleAligment(False)

testCases.compareTestCases()
