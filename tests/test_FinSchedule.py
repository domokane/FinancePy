###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date, set_date_format, DateFormatTypes
from financepy.utils.calendar import CalendarTypes, Calendar
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes


termination_dateAdjust = True


def check_frequency(schedule, start=0):
    dates = schedule._adjusted_dates
    diff_d1d2 = (schedule._termination_date - schedule._effective_date) / 365.0

    for i in range(start, len(dates) - 2):
        diff = (dates[i+1] - dates[i]) / 365.0
        err = diff - (diff_d1d2 / (len(dates) - 1))
        print(err)
        assert round(err, 1) == 0.0


def test_backward_frequencies():
    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD

    freq_type = FrequencyTypes.SEMI_ANNUAL
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 5
    check_frequency(schedule)

    freq_type = FrequencyTypes.QUARTERLY
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 9
    check_frequency(schedule)

    freq_type = FrequencyTypes.MONTHLY
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 25
    check_frequency(schedule)


def test_forward_frequencies():
    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.ANNUAL
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 3
    check_frequency(schedule)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 5
    check_frequency(schedule)

    freq_type = FrequencyTypes.MONTHLY
    schedule = Schedule(d1,
                        d2,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type,
                        termination_dateAdjust)
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 25
    check_frequency(schedule)


def test_backward_front_stub():
    # BACKWARD SHORT STUB AT FRONT
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
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 9
    check_frequency(schedule, start=1)

    # BACKWARD SUPER SHORT STUB AT FRONT
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
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 9
    check_frequency(schedule, start=1)

def test_forward_end_stub():
    # FORWARD SHORT STUB AT END
    termination_dateAdjust = True

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
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 5
    check_frequency(schedule)

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
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 9
    check_frequency(schedule)

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
    adjusted_dates = schedule._adjusted_dates
    assert len(adjusted_dates) == 5
    check_frequency(schedule)
