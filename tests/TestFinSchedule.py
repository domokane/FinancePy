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

def dumpSchedule(desc, schedule):

    testCases.banner("=======================================================")
    testCases.banner(desc)
    testCases.banner("=======================================================")
    testCases.header("OBJ")
    testCases.print(schedule)

    testCases.header("NUM", "TYPE", "DATE", "YEAR", "DIFF")

    numFlows = len(schedule._adjustedDates)
    effDate = schedule._adjustedDates[0]
    years = 0.0
    diff = 0.0
    testCases.print(0, "EFCT DATE", str(effDate), years, diff)
    
    prevDate = schedule._adjustedDates[0]
    for iFlow in range(1, numFlows-1):
        adjustedDate = schedule._adjustedDates[iFlow]
        years = (adjustedDate - effDate) / 365.0
        diff = (adjustedDate - prevDate) / 365.0
        testCases.print(iFlow, "FLOW DATE", str(adjustedDate), years, diff)
        prevDate = adjustedDate

    termDate = schedule._adjustedDates[-1]
    years = (termDate - effDate) / 365.0
    diff = (termDate - prevDate) / 365.0

    testCases.print(numFlows-1, "TERM DATE", str(termDate), years, diff)

############################################################################### 
   
def test_FinSchedule():

    ###########################################################################
    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################
    
    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    terminationDateAdjust = True

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("BACKWARD SEMI-ANNUAL FREQUENCY", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.QUARTERLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("BACKWARD QUARTERLY FREQUENCY", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("BACKWARD MONTHLY FREQUENCY", schedule)

    ###########################################################################
    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)
    
    dumpSchedule("FORWARD ANNUAL", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
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
    
    dumpSchedule("FORWARD SEMI-ANNUAL", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.MONTHLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)
    
    dumpSchedule("FORWARD MONTHLY", schedule)

    ###########################################################################
    # BACKWARD SHORT STUB AT FRONT
    ###########################################################################

    d1 = FinDate(20, 8, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.QUARTERLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)
    dumpSchedule("BACKWARD GEN WITH SHORT END STUB", schedule)

    ###########################################################################
    # BACKWARD SUPER SHORT STUB AT FRONT
    ###########################################################################

    d1 = FinDate(19, 9, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.QUARTERLY
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("BACKWARD GEN WITH VERY SHORT END STUB", schedule)

    ###########################################################################
    # FORWARD SHORT STUB AT END
    ###########################################################################

    d1 = FinDate(20, 8, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("FORWARD GEN WITH END STUB", schedule)

    d1 = FinDate(19, 9 , 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.QUARTERLY
    calendarType = FinCalendarTypes.TARGET
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.FORWARD

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType)

    dumpSchedule("FORWARD GEN WITH VERY SHORT END STUB", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    terminationDateAdjust = True
    
    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust)

    dumpSchedule("TERMINATION DATE ADJUSTED", schedule)

    d1 = FinDate(20, 6, 2018)
    d2 = FinDate(20, 6, 2020)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    terminationDateAdjust = True
    eomFlag = True

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust,
                           eomFlag)

    dumpSchedule("END OF MONTH - NOT EOM TERM DATE - USING MOD FOLL", schedule)

    d1 = FinDate(30, 6, 2018)
    d2 = FinDate(30, 6, 2020)
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    calendarType = FinCalendarTypes.WEEKEND
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    terminationDateAdjust = True
    eomFlag = True

    schedule = FinSchedule(d1,
                           d2,
                           freqType,
                           calendarType,
                           busDayAdjustType,
                           dateGenRuleType,
                           terminationDateAdjust,
                           eomFlag)

    dumpSchedule("END OF MONTH - EOM TERM DATE - USING MOD FOLL", schedule)
    
###############################################################################

def test_FinScheduleAlignment(eomFlag):
        
    valuationDate = FinDate(29, 3, 2005)
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

def test_FinScheduleAlignmentLeapYearEOM():
    ''' Effective date on leap year.'''
    
    valuationDate = FinDate(26, 2, 2006)
    effDate = valuationDate.addTenor("2D")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    calendarType = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")
    eomFlag = True
    
    sched1 = FinSchedule(effDate,
                         matDate1,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)
        
    sched2 = FinSchedule(effDate,
                         matDate2,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)

#    print(sched1._adjustedDates)
#    print(sched2._adjustedDates[:len(sched1._adjustedDates)])

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == eomFlag)

###############################################################################

def test_FinScheduleAlignmentLeapYearNotEOM():
    ''' Effective date on leap year. Not EOM. '''
    
    eomFlag = False

    valuationDate = FinDate(26, 2, 2006)
    effDate = valuationDate.addTenor("2D")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    calendarType = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")

#    print(matDate1, matDate2) 

    sched1 = FinSchedule(effDate,
                         matDate1,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)
        
    sched2 = FinSchedule(effDate,
                         matDate2,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)

#    print(sched1._adjustedDates)
#    print(sched2._adjustedDates[:len(sched1._adjustedDates)])

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == True)

###############################################################################

def test_FinScheduleAlignmentEff31():
    ''' EOM schedule so all unadjusted dates fall on month end.'''
    
    eomFlag = True
    valuationDate = FinDate(29, 7, 2006)
    effDate = valuationDate.addTenor("2D")
    freqType = FinFrequencyTypes.SEMI_ANNUAL
    busDayAdjustType = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    dateGenRuleType = FinDateGenRuleTypes.BACKWARD
    calendarType = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")
    
#    print(matDate1, matDate2)

    sched1 = FinSchedule(effDate,
                         matDate1,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)
        
    sched2 = FinSchedule(effDate,
                         matDate2,
                         freqType,
                         calendarType,
                         busDayAdjustType,
                         dateGenRuleType,
                         adjustTerminationDate, 
                         eomFlag)

#    print(sched1._adjustedDates)
#    print(sched2._adjustedDates[:len(sched1._adjustedDates)])

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == True)

###############################################################################

test_FinSchedule()
test_FinScheduleAlignment(True)
test_FinScheduleAlignment(False)

test_FinScheduleAlignmentLeapYearEOM()
test_FinScheduleAlignmentLeapYearNotEOM()

test_FinScheduleAlignmentEff31()

testCases.compareTestCases()

setDateFormatType(FinDateFormatTypes.UK_LONGEST)
