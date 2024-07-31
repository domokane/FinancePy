###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.schedule import Schedule
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes, Calendar
from financepy.utils.date import Date, set_date_format, DateFormatTypes

test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

set_date_format(DateFormatTypes.UK_LONGEST)


def dumpSchedule(desc, schedule):

    test_cases.banner("=======================================================")
    test_cases.banner(desc)
    test_cases.banner("=======================================================")
    test_cases.header("OBJ")
    test_cases.print(schedule)

    test_cases.header("NUM", "TYPE", "DATE", "YEAR", "DIFF")

    num_flows = len(schedule.adjusted_dts)
    effDate = schedule.adjusted_dts[0]
    years = 0.0
    diff = 0.0
    test_cases.print(0, "EFCT DATE", str(effDate), years, diff)

    prev_dt = schedule.adjusted_dts[0]
    for i_flow in range(1, num_flows-1):
        adjustedDate = schedule.adjusted_dts[i_flow]
        years = (adjustedDate - effDate) / 365.0
        diff = (adjustedDate - prev_dt) / 365.0
        test_cases.print(i_flow, "FLOW DATE", str(adjustedDate), years, diff)
        prev_dt = adjustedDate

    termDate = schedule.adjusted_dts[-1]
    years = (termDate - effDate) / 365.0
    diff = (termDate - prev_dt) / 365.0

    test_cases.print(num_flows-1, "TERM DATE", str(termDate), years, diff)

###############################################################################


def test_FinSchedule():

    ###########################################################################
    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dtAdjust = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("BACKWARD SEMI-ANNUAL FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("BACKWARD QUARTERLY FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("BACKWARD MONTHLY FREQUENCY", schedule)

    ###########################################################################
    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    ###########################################################################

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("FORWARD ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type)

    dumpSchedule("FORWARD SEMI-ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("FORWARD MONTHLY", schedule)

    ###########################################################################
    # BACKWARD SHORT STUB AT FRONT
    ###########################################################################

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)
    dumpSchedule("BACKWARD GEN WITH SHORT END STUB", schedule)

    ###########################################################################
    # BACKWARD SUPER SHORT STUB AT FRONT
    ###########################################################################

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("BACKWARD GEN WITH VERY SHORT END STUB", schedule)

    ###########################################################################
    # FORWARD SHORT STUB AT END
    ###########################################################################

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("FORWARD GEN WITH END STUB", schedule)

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type)

    dumpSchedule("FORWARD GEN WITH VERY SHORT END STUB", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dtAdjust = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust)

    dumpSchedule("TERMINATION DATE ADJUSTED", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dtAdjust = True
    eomFlag = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust,
                        eomFlag)

    dumpSchedule("END OF MONTH - NOT EOM TERM DATE - USING MOD FOLL", schedule)

    d1 = Date(30, 6, 2018)
    d2 = Date(30, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dtAdjust = True
    eomFlag = True

    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        cal_type,
                        bd_type,
                        dg_type,
                        termination_dtAdjust,
                        eomFlag)

    dumpSchedule("END OF MONTH - EOM TERM DATE - USING MOD FOLL", schedule)

    # PROBLEM WITH THIS ONE AS DATES COLLIDE BUT REMOVE FIRST ONE
    schedule = Schedule(Date(28, 4, 2023),
                        Date(30, 4, 2024),
                        FrequencyTypes.ANNUAL,
                        CalendarTypes.UNITED_STATES,
                        BusDayAdjustTypes.MODIFIED_FOLLOWING,
                        DateGenRuleTypes.BACKWARD)

#    print(schedule)
#    print(schedule.adjusted_dts)

###############################################################################


def test_FinScheduleAlignment(eomFlag):

    value_dt = Date(29, 3, 2005)
    effDate = value_dt.add_tenor("2d")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = False

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

#    print(matDate1)
#    print(matDate2)

    myCal = Calendar(cal_type)

    adjustedMatDate1 = myCal.adjust(matDate1, bd_type)
    adjustedMatDate2 = myCal.adjust(matDate2, bd_type)

#    print(adjustedMatDate1)
#    print(adjustedMatDate2)

    sched1 = Schedule(effDate,
                      adjustedMatDate1,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

#    print(sched1)

    sched2 = Schedule(effDate,
                      adjustedMatDate2,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

#    print(sched1.adjusted_dts[-1])
#    print(sched2.adjusted_dts[len(sched1.adjusted_dts)-1])

# THIS TEST IS NO LONGER CORRECT AS I HAVE CHANGED THE  LOGIC TO STEP IN MULTIPLES

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts)-1])


#    print(compare, eomFlag)
#    assert(compare == eomFlag)

###############################################################################


def test_FinScheduleAlignmentLeapYearEOM():
    """ Effective date on leap year."""

    value_dt = Date(26, 2, 2006)
    effDate = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")
    eomFlag = True

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

#    print(sched1.adjusted_dts)
#    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts)-1])
    assert(compare == eomFlag)

###############################################################################


def test_FinScheduleAlignmentLeapYearNotEOM():
    """ Effective date on leap year. Not EOM. """

    eomFlag = False

    value_dt = Date(26, 2, 2006)
    effDate = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

#    print(matDate1, matDate2)

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

#    print(sched1.adjusted_dts)
#    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts)-1])
    assert( compare == True )

###############################################################################


def test_FinScheduleAlignmentEff31():
    """ EOM schedule so all unadjusted dates fall on month end."""

    eomFlag = True
    value_dt = Date(29, 7, 2006)
    effDate = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    matDate1 = effDate.add_tenor("4Y")
    matDate2 = effDate.add_tenor("50Y")

#    print(matDate1, matDate2)

    sched1 = Schedule(effDate,
                      matDate1,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

    sched2 = Schedule(effDate,
                      matDate2,
                      freq_type,
                      cal_type,
                      bd_type,
                      dg_type,
                      adjust_termination_dt,
                      eomFlag)

#    print(sched1.adjusted_dts)
#    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts)-1])
    assert(compare == True)

###############################################################################


test_FinSchedule()
test_FinScheduleAlignment(True)
test_FinScheduleAlignment(False)

test_FinScheduleAlignmentLeapYearEOM()
test_FinScheduleAlignmentLeapYearNotEOM()

test_FinScheduleAlignmentEff31()

test_cases.compareTestCases()

set_date_format(DateFormatTypes.UK_LONGEST)
