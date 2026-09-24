# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes


termination_date_adjust = True

########################################################################################


def check_frequency(schedule, start=0):

    dates = schedule.adjusted_dts
    diff_d1d2 = (schedule.termination_dt - schedule.effective_dt) / 365.0

    for i in range(start, len(dates) - 2):
        diff = (dates[i + 1] - dates[i]) / 365.0
        err = diff - (diff_d1d2 / (len(dates) - 1))
        print(err)
        assert round(err, 1) == 0.0


########################################################################################


def test_backward_frequencies():

    # BACKWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    freq_type = FrequencyTypes.SEMI_ANNUAL
    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 5
    check_frequency(schedule)

    freq_type = FrequencyTypes.QUARTERLY
    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 9
    check_frequency(schedule)

    freq_type = FrequencyTypes.MONTHLY
    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 25
    check_frequency(schedule)


########################################################################################


def test_forward_frequencies():

    # FORWARD SCHEDULES TESTING DIFFERENT FREQUENCIES
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.ANNUAL
    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 3
    check_frequency(schedule)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    schedule = Schedule(d1, d2, freq_type, cal_type, bd_type, dg_type)
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 5
    check_frequency(schedule)

    freq_type = FrequencyTypes.MONTHLY
    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 25
    check_frequency(schedule)


########################################################################################


def test_backward_front_stub():

    # BACKWARD SHORT STUB AT FRONT
    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 9
    check_frequency(schedule, start=1)

    # BACKWARD SUPER SHORT STUB AT FRONT
    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 9
    check_frequency(schedule, start=1)


########################################################################################


def test_forward_end_stub():

    # FORWARD SHORT STUB AT END
    termination_date_adjust = True

    d1 = Date(20, 8, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 5
    check_frequency(schedule)

    d1 = Date(19, 9, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.QUARTERLY
    cal_type = CalendarTypes.TARGET
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    schedule = Schedule(d1, d2, freq_type, cal_type, bd_type, dg_type)
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 9
    check_frequency(schedule)

    d1 = Date(20, 6, 2018)
    d2 = Date(20, 6, 2020)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD
    termination_date_adjust = True

    schedule = Schedule(
        d1, d2, freq_type, cal_type, bd_type, dg_type, termination_date_adjust
    )
    adjusted_dts = schedule.adjusted_dts
    assert len(adjusted_dts) == 5
    check_frequency(schedule)


########################################################################################


def test_generate_is_idempotent_when_termination_date_is_adjusted():
    """Calling generate() again must give the same dates as the first call, which
    the constructor makes. The termination date falls on a Saturday so that it is
    business day adjusted; stepping back from the adjusted date would shift the
    unadjusted dates and add a spurious short stub."""
    d1 = Date(1, 7, 2026)
    d2 = Date(1, 7, 2028)  # Saturday
    schedule = Schedule(
        d1,
        d2,
        FrequencyTypes.SEMI_ANNUAL,
        CalendarTypes.WEEKEND,
        BusDayAdjustTypes.FOLLOWING,
        DateGenRuleTypes.BACKWARD,
        termination_date_adjust,
    )
    first = list(schedule.adjusted_dts)
    expected = [
        Date(1, 7, 2026),
        Date(1, 1, 2027),
        Date(1, 7, 2027),
        Date(3, 1, 2028),
        Date(3, 7, 2028),
    ]
    assert first == expected
    assert schedule.generate() == expected
    assert schedule.generate() == expected
    assert schedule.schedule_dts() == expected
    assert schedule.termination_dt == Date(3, 7, 2028)

    # Products build a schedule and call generate() on it, e.g. the swap legs
    schedule = Schedule(d1, d2, FrequencyTypes.QUARTERLY)
    assert schedule.generate() == list(schedule.adjusted_dts)
