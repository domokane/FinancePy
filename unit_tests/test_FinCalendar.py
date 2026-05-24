# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.date import Date

bus_days_in_decade = {
    CalendarTypes.NONE: 3653,
    CalendarTypes.WEEKEND: 2609,
    CalendarTypes.AUSTRALIA: 2517,
    CalendarTypes.CANADA: 2502,
    CalendarTypes.FRANCE: 2507,
    CalendarTypes.GERMANY: 2529,
    CalendarTypes.ITALY: 2519,
    CalendarTypes.JAPAN: 2466,
    CalendarTypes.NEW_ZEALAND: 2520,
    CalendarTypes.NORWAY: 2526,
    CalendarTypes.SWEDEN: 2514,
    CalendarTypes.SWITZERLAND: 2530,
    CalendarTypes.TARGET: 2562,
    CalendarTypes.UNITED_STATES: 2507,
    CalendarTypes.UNITED_KINGDOM: 2527,
}

MONTHS = {
    "JAN": 1,
    "FEB": 2,
    "MAR": 3,
    "APR": 4,
    "MAY": 5,
    "JUN": 6,
    "JUL": 7,
    "AUG": 8,
    "SEP": 9,
    "OCT": 10,
    "NOV": 11,
    "DEC": 12,
}


def date_from_string(s):
    _, d, m, y = s.split()
    return Date(int(d), MONTHS[m], int(y))


########################################################################################


import pytest


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_round_trip(cal_type):

    cal = Calendar(cal_type)

    start_dates = [
        Date(3, 1, 2020),
        Date(28, 2, 2020),
        Date(24, 12, 2021),
        Date(31, 12, 2024),
    ]

    steps = [0, 1, 2, 5, 10, 50, 250, -1, -2, -5, -10, -50, -250]

    for start in start_dates:

        if not cal.is_business_day(start):
            continue

        for n in steps:
            end = cal.add_business_days(start, n)
            back = cal.add_business_days(end, -n)

            assert back == start, (
                f"{cal_type}: round-trip failed: "
                f"{start} + {n} business days = {end}, "
                f"then back {-n} gives {back}"
            )


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_from_non_business_day_lands_consistently(cal_type):

    cal = Calendar(cal_type)
    start = Date(24, 12, 2021)

    if cal.is_business_day(start):
        return

    next_bd = cal.add_business_days(start, 1)
    prev_bd = cal.add_business_days(start, -1)

    assert cal.is_business_day(next_bd)
    assert cal.is_business_day(prev_bd)
    assert next_bd > start
    assert prev_bd < start


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_matches_day_by_day_count(cal_type):

    cal = Calendar(cal_type)

    start = Date(3, 1, 2020)

    for n in [1, 2, 5, 10, 50, 100, 250]:

        d = start
        count = 0

        while count < n:
            d = d.add_days(1)
            if cal.is_business_day(d):
                count += 1

        assert cal.add_business_days(start, n) == d, (
            f"{cal_type}: add_business_days disagrees with "
            f"day-by-day business-day counting for n={n}"
        )


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_round_trip_from_business_days(cal_type):

    cal = Calendar(cal_type)

    start_dates = [
        Date(3, 1, 2020),
        Date(28, 2, 2020),
        Date(23, 12, 2021),
        Date(31, 12, 2024),
    ]

    steps = [0, 1, 2, 5, 10, 50, 250, -1, -2, -5, -10, -50, -250]

    for start in start_dates:

        if not cal.is_business_day(start):
            continue

        for n in steps:
            end = cal.add_business_days(start, n)
            back = cal.add_business_days(end, -n)

            assert back == start, (
                f"{cal_type}: {start} + {n} business days = {end}, "
                f"then back gives {back}"
            )


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_round_trip_from_business_days(cal_type):

    cal = Calendar(cal_type)

    start_dates = [
        Date(3, 1, 2020),
        Date(28, 2, 2020),
        Date(23, 12, 2021),
        Date(31, 12, 2024),
    ]

    steps = [0, 1, 2, 5, 10, 50, 250, -1, -2, -5, -10, -50, -250]

    for start in start_dates:

        if not cal.is_business_day(start):
            continue

        for n in steps:
            end = cal.add_business_days(start, n)
            back = cal.add_business_days(end, -n)

            assert back == start, (
                f"{cal_type}: {start} + {n} business days = {end}, "
                f"then back gives {back}"
            )


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_adjust_none_returns_same_date(cal_type):
    cal = Calendar(cal_type)
    dt = Date(4, 7, 2020)

    assert cal.adjust(dt, BusDayAdjustTypes.NONE) == dt


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
@pytest.mark.parametrize("bd_type", list(BusDayAdjustTypes))
def test_adjusted_date_is_business_day_unless_no_adjustment(cal_type, bd_type):
    cal = Calendar(cal_type)
    dt = Date(4, 7, 2020)

    adjusted = cal.adjust(dt, bd_type)

    if cal_type is CalendarTypes.NONE or bd_type is BusDayAdjustTypes.NONE:
        assert adjusted == dt
    else:
        assert cal.is_business_day(adjusted)


def test_us_observed_independence_day_2010():
    cal = Calendar(CalendarTypes.UNITED_STATES)

    assert cal.is_holiday(Date(4, 7, 2010))
    assert cal.is_holiday(Date(5, 7, 2010))
    assert not cal.is_business_day(Date(5, 7, 2010))
    assert cal.adjust(Date(4, 7, 2010), BusDayAdjustTypes.FOLLOWING) == Date(6, 7, 2010)


def test_modified_following_does_not_cross_month():
    cal = Calendar(CalendarTypes.WEEKEND)

    dt = Date(31, 7, 2021)  # Saturday
    adjusted = cal.adjust(dt, BusDayAdjustTypes.MODIFIED_FOLLOWING)

    assert adjusted == Date(30, 7, 2021)


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_zero_business_days_returns_same_date(cal_type):
    cal = Calendar(cal_type)
    dt = Date(24, 12, 2021)

    assert cal.add_business_days(dt, 0) == dt


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_add_business_days_old_matches_new(cal_type):
    cal = Calendar(cal_type)
    start = Date(3, 1, 2020)

    for n in [-20, -5, -1, 0, 1, 5, 20]:
        assert cal.add_business_days(start, n) == cal.add_business_days_old(start, n)


@pytest.mark.parametrize("cal_type", list(CalendarTypes))
def test_get_holiday_list_contains_only_holidays(cal_type):

    cal = Calendar(cal_type)

    for holiday in cal.get_holiday_list(2021):
        dt = date_from_string(holiday)
        assert cal.is_holiday(dt)
