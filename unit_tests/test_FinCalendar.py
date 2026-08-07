# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.date import Date
from financepy.utils.date_format import set_date_format
from financepy.utils.date_format import DateFormatTypes
from financepy.utils.error import FinError

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
    CalendarTypes.LONDON: 2527,
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
    set_date_format(DateFormatTypes.UK_LONGEST)
    for holiday in cal.get_holiday_list(2021):
        dt = date_from_string(holiday)
        assert cal.is_holiday(dt)

@pytest.mark.parametrize("cal_type", list(CalendarTypes))
@pytest.mark.parametrize("bd_type", list(BusDayAdjustTypes))
def test_fast_adjust_matches_adjust(cal_type, bd_type):

    cal = Calendar(cal_type)

    dates = [
        Date(31,1,2021),
        Date(28,2,2021),
        Date(24,12,2021),
        Date(25,12,2021),
        Date(31,12,2021),
    ]

    for dt in dates:
        assert cal.fast_adjust(dt, bd_type) == cal.adjust(dt, bd_type)


########################################################################################
# Joint calendars — a list of calendar types unions their holidays
########################################################################################


JOINT_US_UK = [CalendarTypes.UNITED_STATES, CalendarTypes.UNITED_KINGDOM]


def test_joint_calendar_holiday_is_union():

    us = Calendar(CalendarTypes.UNITED_STATES)
    uk = Calendar(CalendarTypes.UNITED_KINGDOM)
    joint = Calendar(JOINT_US_UK)

    july_4_2023 = Date(4, 7, 2023)  # US Independence Day, UK business day
    assert us.is_holiday(july_4_2023)
    assert not uk.is_holiday(july_4_2023)
    assert joint.is_holiday(july_4_2023)

    aug_28_2023 = Date(28, 8, 2023)  # UK late summer bank holiday, US business day
    assert not us.is_holiday(aug_28_2023)
    assert uk.is_holiday(aug_28_2023)
    assert joint.is_holiday(aug_28_2023)

    plain_day = Date(6, 6, 2023)
    assert not joint.is_holiday(plain_day)
    assert joint.is_business_day(plain_day)


def test_joint_business_day_is_intersection_over_full_years():

    us = Calendar(CalendarTypes.UNITED_STATES)
    uk = Calendar(CalendarTypes.UNITED_KINGDOM)
    joint = Calendar(JOINT_US_UK)

    dt = Date(1, 1, 2020)
    end = Date(1, 1, 2024)

    while dt < end:
        assert joint.is_business_day(dt) == (
            us.is_business_day(dt) and uk.is_business_day(dt)
        ), f"joint disagrees with single-calendar intersection on {dt}"
        dt = dt.add_days(1)


def test_joint_adjust_skips_holidays_of_both_markets():

    # Good Friday 2021-04-02 and Easter Monday 2021-04-05 are UK bank
    # holidays but US business days, so a joint US+UK calendar must jump
    # from the Friday all the way to Tuesday 2021-04-06
    us = Calendar(CalendarTypes.UNITED_STATES)
    joint = Calendar(JOINT_US_UK)

    good_friday = Date(2, 4, 2021)
    assert us.adjust(good_friday, BusDayAdjustTypes.FOLLOWING) == good_friday
    assert joint.adjust(good_friday, BusDayAdjustTypes.FOLLOWING) == Date(6, 4, 2021)


@pytest.mark.parametrize("bd_type", list(BusDayAdjustTypes))
def test_joint_fast_adjust_matches_adjust(bd_type):

    joint = Calendar(JOINT_US_UK)

    dates = [
        Date(2, 4, 2021),
        Date(31, 1, 2021),
        Date(24, 12, 2021),
        Date(25, 12, 2021),
        Date(31, 12, 2021),
    ]

    for dt in dates:
        assert joint.fast_adjust(dt, bd_type) == joint.adjust(dt, bd_type)


def test_joint_add_business_days_round_trip():

    joint = Calendar(JOINT_US_UK)
    start = Date(3, 1, 2020)

    for n in [0, 1, 5, 50, 250, -1, -5, -50, -250]:
        end = joint.add_business_days(start, n)
        back = joint.add_business_days(end, -n)
        assert back == start


def test_joint_holiday_list_is_union_of_lists():

    us = Calendar(CalendarTypes.UNITED_STATES)
    uk = Calendar(CalendarTypes.UNITED_KINGDOM)
    joint = Calendar(JOINT_US_UK)

    year = 2023
    us_holidays = set(us.get_holiday_list(year))
    uk_holidays = set(uk.get_holiday_list(year))

    assert set(joint.get_holiday_list(year)) == us_holidays | uk_holidays


def test_joint_single_entry_and_duplicates_collapse():

    single = Calendar([CalendarTypes.UNITED_STATES])
    assert single.cal_type == CalendarTypes.UNITED_STATES

    duplicate = Calendar([CalendarTypes.UNITED_STATES, CalendarTypes.UNITED_STATES])
    assert duplicate.cal_type == CalendarTypes.UNITED_STATES


def test_joint_str_and_repr():

    joint = Calendar(JOINT_US_UK)
    assert str(joint) == "JOINT(UNITED_STATES,UNITED_KINGDOM)"
    assert repr(joint) == str(joint)


def test_joint_invalid_inputs_raise():

    with pytest.raises(FinError):
        Calendar([])

    with pytest.raises(FinError):
        Calendar([CalendarTypes.UNITED_STATES, "TARGET"])
