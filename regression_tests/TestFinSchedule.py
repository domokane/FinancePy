# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.schedule import Schedule
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.calendar import CalendarTypes, Calendar
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date

from FinTestCases import FinTestCases
from FinTestCases import global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

set_date_format(DateFormatTypes.UK_LONGEST)

########################################################################################


def dump_schedule(desc, schedule):

    test_cases.banner("=======================================================")
    test_cases.banner(desc)
    test_cases.banner("=======================================================")
    test_cases.header("OBJ")
    test_cases.print(schedule)

    test_cases.header("NUM", "TYPE", "DATE", "YEAR", "DIFF")

    num_flows = len(schedule.adjusted_dts)
    eff_date = schedule.adjusted_dts[0]
    years = 0.0
    diff = 0.0
    test_cases.print(0, "EFCT DATE", str(eff_date), years, diff)

    prev_dt = schedule.adjusted_dts[0]
    for i_flow in range(1, num_flows - 1):
        adjusted_date = schedule.adjusted_dts[i_flow]
        years = (adjusted_date - eff_date) / 365.0
        diff = (adjusted_date - prev_dt) / 365.0
        test_cases.print(i_flow, "FLOW DATE", str(adjusted_date), years, diff)
        prev_dt = adjusted_date

    term_date = schedule.adjusted_dts[-1]
    years = (term_date - eff_date) / 365.0
    diff = (term_date - prev_dt) / 365.0

    test_cases.print(num_flows - 1, "TERM DATE", str(term_date), years, diff)


########################################################################################


def test_fin_schedule():

    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dt_adjust = True

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("BACKWARD SEMI-ANNUAL FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("BACKWARD QUARTERLY FREQUENCY", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("BACKWARD MONTHLY FREQUENCY", schedule)

    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("FORWARD ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(d1, d2, freq_type, cal_type, bd_type, dg_type)

    dump_schedule("FORWARD SEMI-ANNUAL", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.MONTHLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("FORWARD MONTHLY", schedule)

    # BACKWARD SHORT STUB AT FRONT

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )
    dump_schedule("BACKWARD GEN WITH SHORT END STUB", schedule)

    # BACKWARD SUPER SHORT STUB AT FRONT

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("BACKWARD GEN WITH VERY SHORT END STUB", schedule)

    # FORWARD SHORT STUB AT END

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("FORWARD GEN WITH END STUB", schedule)

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1, d2, freq_type, cal_type, bd_type, dg_type)

    dump_schedule("FORWARD GEN WITH VERY SHORT END STUB", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dt_adjust = True

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_dt_adjust
    )

    dump_schedule("TERMINATION DATE ADJUSTED", schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dt_adjust = True
    eom_flag = True

    schedule = Schedule(
        d1,
        d2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        termination_dt_adjust,
        eom_flag,
    )

    dump_schedule("END OF MONTH - NOT EOM TERM DATE - USING MOD FOLL", schedule)

    d1 = Date(30, 6, 2018)
    d2 = Date(30, 6, 2020)
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_dt_adjust = True
    eom_flag = True

    schedule = Schedule(
        d1,
        d2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        termination_dt_adjust,
        eom_flag,
    )

    dump_schedule("END OF MONTH - EOM TERM DATE - USING MOD FOLL", schedule)

    # PROBLEM WITH THIS ONE AS DATES COLLIDE BUT REMOVE FIRST ONE
    schedule = Schedule(
        Date(28, 4, 2023),
        Date(30, 4, 2024),
        FrequencyTypes.ANNUAL,
        CalendarTypes.UNITED_STATES,
        BusDayAdjustTypes.MODIFIED_FOLLOWING,
        DateGenRuleTypes.BACKWARD,
    )


#    print(schedule)
#    print(schedule.adjusted_dts)

########################################################################################


def test_fin_schedule_alignment(eom_flag):

    value_dt = Date(29, 3, 2005)
    eff_date = value_dt.add_tenor("2d")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = False

    mat_date1 = eff_date.add_tenor("4Y")
    mat_date2 = eff_date.add_tenor("50Y")

    #    print(mat_date1)
    #    print(mat_date2)

    my_cal = Calendar(cal_type)

    adjusted_mat_date1 = my_cal.adjust(mat_date1, bd_type)
    adjusted_mat_date2 = my_cal.adjust(mat_date2, bd_type)

    #    print(adjusted_mat_date1)
    #    print(adjusted_mat_date2)

    sched1 = Schedule(
        eff_date,
        adjusted_mat_date1,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    #    print(sched1)

    sched2 = Schedule(
        eff_date,
        adjusted_mat_date2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    #    print(sched1.adjusted_dts[-1])
    #    print(sched2.adjusted_dts[len(sched1.adjusted_dts)-1])

    # THIS TEST IS NO LONGER CORRECT AS I HAVE CHANGED THE  LOGIC TO STEP IN MULTIPLES

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts) - 1]
    )


#    print(compare, eom_flag)
#    assert(compare == eom_flag)

########################################################################################


def test_fin_schedule_alignment_leap_year_eom():
    """Effective date on leap year."""

    value_dt = Date(26, 2, 2006)
    eff_date = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    mat_date1 = eff_date.add_tenor("4Y")
    mat_date2 = eff_date.add_tenor("50Y")
    eom_flag = True

    sched1 = Schedule(
        eff_date,
        mat_date1,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    sched2 = Schedule(
        eff_date,
        mat_date2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    #    print(sched1.adjusted_dts)
    #    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts) - 1]
    )
    assert compare == eom_flag


########################################################################################


def test_fin_schedule_alignment_leap_year_not_eom():
    """Effective date on leap year. Not EOM."""

    eom_flag = False

    value_dt = Date(26, 2, 2006)
    eff_date = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    mat_date1 = eff_date.add_tenor("4Y")
    mat_date2 = eff_date.add_tenor("50Y")

    #    print(mat_date1, mat_date2)

    sched1 = Schedule(
        eff_date,
        mat_date1,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    sched2 = Schedule(
        eff_date,
        mat_date2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    #    print(sched1.adjusted_dts)
    #    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts) - 1]
    )
    assert compare is True


########################################################################################


def test_fin_schedule_alignment_eff31():
    """EOM schedule so all unadjusted dates fall on month end."""

    eom_flag = True
    value_dt = Date(29, 7, 2006)
    eff_date = value_dt.add_tenor("2D")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    cal_type = CalendarTypes.UNITED_STATES
    adjust_termination_dt = True

    mat_date1 = eff_date.add_tenor("4Y")
    mat_date2 = eff_date.add_tenor("50Y")

    #    print(mat_date1, mat_date2)

    sched1 = Schedule(
        eff_date,
        mat_date1,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    sched2 = Schedule(
        eff_date,
        mat_date2,
        freq_type,
        cal_type,
        bd_type,
        dg_type,
        adjust_termination_dt,
        eom_flag,
    )

    #    print(sched1.adjusted_dts)
    #    print(sched2.adjusted_dts[:len(sched1.adjusted_dts)])

    compare = (
        sched1.adjusted_dts[-1] == sched2.adjusted_dts[len(sched1.adjusted_dts) - 1]
    )
    assert compare is True


########################################################################################

test_fin_schedule()
test_fin_schedule_alignment(True)
test_fin_schedule_alignment(False)

test_fin_schedule_alignment_leap_year_eom()
test_fin_schedule_alignment_leap_year_not_eom()

test_fin_schedule_alignment_eff31()

test_cases.compare_test_cases()

set_date_format(DateFormatTypes.UK_LONGEST)
