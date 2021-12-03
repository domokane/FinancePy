###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.schedule import Schedule
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes, Calendar
from financepy.utils.date import Date, set_date_format, DateFormatTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

set_date_format(DateFormatTypes.UK_LONGEST)


def dumpSchedule(desc, schedule):

    testCases.banner("=======================================================")
    testCases.banner(desc)
    testCases.banner("=======================================================")
    testCases.header("OBJ")
    testCases.print(schedule)

    testCases.header("NUM", "TYPE", "DATE", "YEAR", "DIFF")

    num_flows = len(schedule._adjusted_dates)
    effDate = schedule._adjusted_dates[0]
    years = 0.0
    diff = 0.0
    testCases.print(0, "EFCT DATE", str(effDate), years, diff)

    prev_date = schedule._adjusted_dates[0]
    for iFlow in range(1, num_flows-1):
        adjustedDate = schedule._adjusted_dates[iFlow]
        years = (adjustedDate - effDate) / 365.0
        diff = (adjustedDate - prev_date) / 365.0
        testCases.print(iFlow, "FLOW DATE", str(adjustedDate), years, diff)
        prev_date = adjustedDate

    termDate = schedule._adjusted_dates[-1]
    years = (termDate - effDate) / 365.0
    diff = (termDate - prev_date) / 365.0

    testCases.print(num_flows-1, "TERM DATE", str(termDate), years, diff)

###############################################################################


def test_FinSchedule():

    ###########################################################################
    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
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
    freq_type = FrequencyTypes.QUARTERLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

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
    freq_type = FrequencyTypes.MONTHLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

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
    freq_type = FrequencyTypes.ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD

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
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    dumpSchedule("FORWARD SEMI-ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.MONTHLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

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
    freq_type = FrequencyTypes.QUARTERLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

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
    freq_type = FrequencyTypes.QUARTERLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

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
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD

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
    freq_type = FrequencyTypes.QUARTERLY
    calendar_type = CalendarTypes.TARGET
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    dumpSchedule("FORWARD GEN WITH VERY SHORT END STUB", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
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
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
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
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
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
    effDate = valuation_date.add_tenor("2d")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    calendar_type = CalendarTypes.UNITED_STATES
    adjust_termination_date = False

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

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
                      adjust_termination_date,
                      eomFlag)

#    print(sched1)

    sched2 = Schedule(effDate,
                      adjustedMatDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

#    print(sched1._adjusted_dates[-1])
#    print(sched2._adjusted_dates[len(sched1._adjusted_dates)-1])

# THIS TEST IS NO LONGER CORRECT AS I HAVE CHANGED THE  LOGIC TO STEP IN MULTIPLES

    compare = (
        sched1._adjusted_dates[-1] == sched2._adjusted_dates[len(sched1._adjusted_dates)-1])


#    print(compare, eomFlag)
#    assert(compare == eomFlag)

###############################################################################


def test_FinScheduleAlignmentLeapYearEOM():
    """ Effective date on leap year."""

    valuation_date = Date(26, 2, 2006)
    effDate = valuation_date.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    calendar_type = CalendarTypes.UNITED_STATES
    adjust_termination_date = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")
    eomFlag = True

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

#    print(sched1._adjusted_dates)
#    print(sched2._adjusted_dates[:len(sched1._adjusted_dates)])

    compare = (
        sched1._adjusted_dates[-1] == sched2._adjusted_dates[len(sched1._adjusted_dates)-1])
    assert(compare == eomFlag)

###############################################################################


def test_FinScheduleAlignmentLeapYearNotEOM():
    """ Effective date on leap year. Not EOM. """

    eomFlag = False

    valuation_date = Date(26, 2, 2006)
    effDate = valuation_date.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    calendar_type = CalendarTypes.UNITED_STATES
    adjust_termination_date = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

#    print(matDate1, matDate2)

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

#    print(sched1._adjusted_dates)
#    print(sched2._adjusted_dates[:len(sched1._adjusted_dates)])

    compare = (
        sched1._adjusted_dates[-1] == sched2._adjusted_dates[len(sched1._adjusted_dates)-1])
    assert(compare == True)

###############################################################################


def test_FinScheduleAlignmentEff31():
    """ EOM schedule so all unadjusted dates fall on month end."""

    eomFlag = True
    valuation_date = Date(29, 7, 2006)
    effDate = valuation_date.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    calendar_type = CalendarTypes.UNITED_STATES
    adjust_termination_date = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

#    print(matDate1, matDate2)

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      calendar_type,
                      bus_day_adjust_type,
                      date_gen_rule_type,
                      adjust_termination_date,
                      eomFlag)

#    print(sched1._adjusted_dates)
#    print(sched2._adjusted_dates[:len(sched1._adjusted_dates)])

    compare = (
        sched1._adjusted_dates[-1] == sched2._adjusted_dates[len(sched1._adjusted_dates)-1])
    assert(compare == True)

###############################################################################


test_FinSchedule()
test_FinScheduleAlignment(True)
test_FinScheduleAlignment(False)

test_FinScheduleAlignmentLeapYearEOM()
test_FinScheduleAlignmentLeapYearNotEOM()

test_FinScheduleAlignmentEff31()

testCases.compareTestCases()

set_date_format(DateFormatTypes.UK_LONGEST)
