###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Calendar import FinBusDayAdjustTypes
from financepy.utils.Calendar import FinDateGenRuleTypes
from financepy.utils.Schedule import Schedule
from financepy.utils.Frequency import FinFrequencyTypes
from financepy.utils.Calendar import FinCalendarTypes, Calendar
from financepy.utils.Date import Date, setDateFormatType, FinDateFormatTypes

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
    
    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    termination_dateAdjust = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("BACKWARD SEMI-ANNUAL FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.QUARTERLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("BACKWARD QUARTERLY FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.MONTHLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("BACKWARD MONTHLY FREQUENCY", schedule)

    ###########################################################################
    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    
    dumpSchedule("FORWARD ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)
    
    dumpSchedule("FORWARD SEMI-ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.MONTHLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    
    dumpSchedule("FORWARD MONTHLY", schedule)

    ###########################################################################
    # BACKWARD SHORT STUB AT FRONT
    ###########################################################################

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.QUARTERLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    dumpSchedule("BACKWARD GEN WITH SHORT END STUB", schedule)

    ###########################################################################
    # BACKWARD SUPER SHORT STUB AT FRONT
    ###########################################################################

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.QUARTERLY
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("BACKWARD GEN WITH VERY SHORT END STUB", schedule)

    ###########################################################################
    # FORWARD SHORT STUB AT END
    ###########################################################################

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("FORWARD GEN WITH END STUB", schedule)

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.QUARTERLY
    calendar_type = FinCalendarTypes.TARGET
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    dumpSchedule("FORWARD GEN WITH VERY SHORT END STUB", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    termination_dateAdjust = True
    
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)

    dumpSchedule("TERMINATION DATE ADJUSTED", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    termination_dateAdjust = True
    eomFlag = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust,
                        eomFlag)

    dumpSchedule("END OF MONTH - NOT EOM TERM DATE - USING MOD FOLL", schedule)

    d1 = Date(30, 6, 2018)
    d2 = Date(30, 6, 2020)
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    calendar_type = FinCalendarTypes.WEEKEND
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    termination_dateAdjust = True
    eomFlag = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust,
                        eomFlag)

    dumpSchedule("END OF MONTH - EOM TERM DATE - USING MOD FOLL", schedule)
    
###############################################################################

def test_FinScheduleAlignment(eomFlag):
        
    valuation_date = Date(29, 3, 2005)
    effDate = valuation_date.addTenor("2d")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    calendar_type = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = False

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")

#    print(matDate1)
#    print(matDate2)

    myCal = Calendar(calendar_type)

    adjustedMatDate1 = myCal.adjust(matDate1, bus_day_adjust_type)
    adjustedMatDate2 = myCal.adjust(matDate2, bus_day_adjust_type)

#    print(adjustedMatDate1)
#    print(adjustedMatDate2)
    
    sched1 = Schedule(effDate,
                      adjustedMatDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)
    
#    print(sched1)
    
    sched2 = Schedule(effDate,
                      adjustedMatDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == eomFlag)

###############################################################################

def test_FinScheduleAlignmentLeapYearEOM():
    """ Effective date on leap year."""
    
    valuation_date = Date(26, 2, 2006)
    effDate = valuation_date.addTenor("2D")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    calendar_type = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")
    eomFlag = True
    
    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)
        
    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)

#    print(sched1._adjustedDates)
#    print(sched2._adjustedDates[:len(sched1._adjustedDates)])

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == eomFlag)

###############################################################################

def test_FinScheduleAlignmentLeapYearNotEOM():
    """ Effective date on leap year. Not EOM. """
    
    eomFlag = False

    valuation_date = Date(26, 2, 2006)
    effDate = valuation_date.addTenor("2D")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    calendar_type = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")

#    print(matDate1, matDate2) 

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)
        
    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)

#    print(sched1._adjustedDates)
#    print(sched2._adjustedDates[:len(sched1._adjustedDates)])

    compare = (sched1._adjustedDates[-1] == sched2._adjustedDates[len(sched1._adjustedDates)-1])
    assert(compare == True)

###############################################################################

def test_FinScheduleAlignmentEff31():
    """ EOM schedule so all unadjusted dates fall on month end."""
    
    eomFlag = True
    valuation_date = Date(29, 7, 2006)
    effDate = valuation_date.addTenor("2D")
    freq_type = FinFrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = FinBusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = FinDateGenRuleTypes.BACKWARD
    calendar_type = FinCalendarTypes.UNITED_STATES
    adjustTerminationDate = True

    matDate1 = effDate.addTenor("4Y")
    matDate2 = effDate.addTenor("50Y")
    
#    print(matDate1, matDate2)

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjustTerminationDate,
                      eomFlag)
        
    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
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
