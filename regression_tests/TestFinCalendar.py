# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date
from financepy.utils.calendar import BusDayAdjustTypes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_calendar():

    set_date_format(DateFormatTypes.US_LONGEST)
    end_dt = Date(31, 12, 2030)

    for cal_type in CalendarTypes:

        test_cases.banner("================================")
        test_cases.banner("================================")

        test_cases.header("CALENDAR", "HOLIDAY")
        test_cases.print("STARTING", cal_type)

        cal = Calendar(cal_type)
        next_dt = Date(31, 12, 2020)

        while next_dt < end_dt:
            next_dt = next_dt.add_days(1)

            if next_dt.d == 1 and next_dt.m == 1:
                test_cases.banner("================================")
            #                print("=========================")

            is_holiday_day = cal.is_holiday(next_dt)
            if is_holiday_day is True:
                test_cases.print(cal, next_dt)
    #                print(cal, next_dt)

    set_date_format(DateFormatTypes.US_LONG)


def test_target():
    cal = Calendar(CalendarTypes.TARGET)

    # --------------------------------------------------------------
    # 2026 TARGET holidays
    # --------------------------------------------------------------

    # New Year's Day
    assert not cal.is_business_day(Date(1, 1, 2026))

    # Good Friday
    assert not cal.is_business_day(Date(3, 4, 2026))

    # Easter Monday
    assert not cal.is_business_day(Date(6, 4, 2026))

    # Labour Day
    assert not cal.is_business_day(Date(1, 5, 2026))

    # Christmas Day
    assert not cal.is_business_day(Date(25, 12, 2026))

    # Boxing Day / 26 December
    assert not cal.is_business_day(Date(26, 12, 2026))

    # --------------------------------------------------------------
    # Ordinary TARGET business days
    # --------------------------------------------------------------

    assert cal.is_business_day(Date(2, 1, 2026))
    assert cal.is_business_day(Date(2, 4, 2026))
    assert cal.is_business_day(Date(7, 4, 2026))
    assert cal.is_business_day(Date(30, 4, 2026))
    assert cal.is_business_day(Date(4, 5, 2026))
    assert cal.is_business_day(Date(24, 12, 2026))
    assert cal.is_business_day(Date(28, 12, 2026))

def test_target_not_national_holidays():
    cal = Calendar(CalendarTypes.TARGET)

    # German Unity Day is not a TARGET closing day.
    # 3 Oct 2025 is Friday.
    assert cal.is_business_day(Date(3, 10, 2025))

    # Bastille Day is not a TARGET closing day.
    assert cal.is_business_day(Date(14, 7, 2026))

    # Assumption Day is not a TARGET holiday.
    # Use a year where it is a weekday.
    assert cal.is_business_day(Date(15, 8, 2025))

    # Armistice Day is not a TARGET holiday.
    assert cal.is_business_day(Date(11, 11, 2026))

def test_fast_adjust():
    for cal_type in CalendarTypes:
        cal = Calendar(cal_type)
        for y in range(2000, 2081):  # DO NOT NEED ENTIRE 200 YEAR RANGE
            dt = Date(1, 1, y)
            end_dt = Date(1, 1, y + 1)

            while dt < end_dt:
                for bd_type in BusDayAdjustTypes:
                    a = cal.adjust(dt, bd_type)
                    b = cal.fast_adjust(dt, bd_type)

                    if a != b:
                        print(
                            f"Mismatch: {cal_type}, {bd_type}, {dt}, {a}, {b}"
                        )

                dt = dt.add_days(1)

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


def test_joint_fast_adjust_matches_adjust():

    joint = Calendar(JOINT_US_UK)

    dates = [
        Date(2, 4, 2021),
        Date(31, 1, 2021),
        Date(24, 12, 2021),
        Date(25, 12, 2021),
        Date(31, 12, 2021),
    ]

    for bd_type in list(BusDayAdjustTypes):
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

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def test_weekend():
    cal = Calendar(CalendarTypes.WEEKEND)

    assert not cal.is_business_day(Date(8, 8, 2026))   # Saturday
    assert not cal.is_business_day(Date(9, 8, 2026))   # Sunday
    assert cal.is_business_day(Date(10, 8, 2026))       # Monday


def test_australia():
    cal = Calendar(CalendarTypes.AUSTRALIA)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(26, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))    # Good Friday
    assert not cal.is_business_day(Date(6, 4, 2026))    # Easter Monday
    assert not cal.is_business_day(Date(25, 4, 2026))   # ANZAC
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_canada():
    cal = Calendar(CalendarTypes.CANADA)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(16, 2, 2026))   # Family Day
    assert not cal.is_business_day(Date(3, 4, 2026))    # Good Friday
    assert not cal.is_business_day(Date(18, 5, 2026))   # Victoria Day
    assert not cal.is_business_day(Date(1, 7, 2026))    # Canada Day
    assert not cal.is_business_day(Date(7, 9, 2026))    # Labour Day
    assert not cal.is_business_day(Date(12, 10, 2026))  # Thanksgiving
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_france():
    cal = Calendar(CalendarTypes.FRANCE)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))    # Easter Monday
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(8, 5, 2026))
    assert not cal.is_business_day(Date(14, 7, 2026))
    assert not cal.is_business_day(Date(11, 11, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_germany():
    cal = Calendar(CalendarTypes.GERMANY)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(14, 5, 2026))   # Ascension
    assert not cal.is_business_day(Date(3, 10, 2025))   # Unity Day
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_italy():
    cal = Calendar(CalendarTypes.ITALY)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(6, 1, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(25, 4, 2026))
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(2, 6, 2026))
    assert not cal.is_business_day(Date(8, 12, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_new_zealand():
    cal = Calendar(CalendarTypes.NEW_ZEALAND)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(2, 1, 2026))
    assert not cal.is_business_day(Date(6, 2, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(25, 4, 2026))
    assert not cal.is_business_day(Date(26, 10, 2026))  # Labour Day
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_norway():
    cal = Calendar(CalendarTypes.NORWAY)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(2, 4, 2026))    # Maundy Thursday
    assert not cal.is_business_day(Date(3, 4, 2026))    # Good Friday
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(14, 5, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_sweden():
    cal = Calendar(CalendarTypes.SWEDEN)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(6, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(6, 6, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_switzerland():
    cal = Calendar(CalendarTypes.SWITZERLAND)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(2, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(14, 5, 2026))
    assert not cal.is_business_day(Date(1, 8, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_united_states():
    cal = Calendar(CalendarTypes.UNITED_STATES)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(19, 1, 2026))
    assert not cal.is_business_day(Date(16, 2, 2026))
    assert not cal.is_business_day(Date(25, 5, 2026))
    assert not cal.is_business_day(Date(19, 6, 2026))
    assert not cal.is_business_day(Date(4, 7, 2026))
    assert not cal.is_business_day(Date(7, 9, 2026))
    assert not cal.is_business_day(Date(26, 11, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_new_york():
    cal = Calendar(CalendarTypes.NEW_YORK)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(19, 1, 2026))
    assert not cal.is_business_day(Date(16, 2, 2026))
    assert not cal.is_business_day(Date(25, 5, 2026))
    assert not cal.is_business_day(Date(19, 6, 2026))
    assert not cal.is_business_day(Date(7, 9, 2026))
    assert not cal.is_business_day(Date(12, 10, 2026))
    assert not cal.is_business_day(Date(11, 11, 2026))
    assert not cal.is_business_day(Date(26, 11, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))


def test_london():
    cal = Calendar(CalendarTypes.LONDON)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(4, 5, 2026))    # Early May
    assert not cal.is_business_day(Date(25, 5, 2026))   # Spring
    assert not cal.is_business_day(Date(31, 8, 2026))   # Summer
    assert not cal.is_business_day(Date(25, 12, 2026))
    assert not cal.is_business_day(Date(28, 12, 2026))  # Boxing observed


def test_united_kingdom():
    london = Calendar(CalendarTypes.LONDON)
    uk = Calendar(CalendarTypes.UNITED_KINGDOM)

    # If these are intentionally equivalent, enforce it.
    dt = Date(1, 1, 2020)
    end = Date(1, 1, 2031)

    while dt < end:
        assert uk.is_business_day(dt) == london.is_business_day(dt), (
            f"UK/London mismatch on {dt}"
        )
        dt = dt.add_days(1)

def test_japan():

    cal = Calendar(CalendarTypes.JAPAN)

    # Jan 2 is not a Japanese statutory holiday
    assert cal.is_business_day(Date(2, 1, 2026))

    # Jan 3, 2026 is Saturday, so not a business day,
    # but it is not a statutory holiday.
    assert not cal.is_business_day(Date(3, 1, 2026))
    assert not cal.is_holiday(Date(3, 1, 2026))

    # Vernal Equinox
    assert not cal.is_business_day(Date(20, 3, 2026))

    # Substitute holiday:
    # Constitution Memorial Day is Sunday May 3,
    # so May 6 becomes substitute holiday after May 4/5.
    assert not cal.is_business_day(Date(6, 5, 2026))

    # Marine Day
    assert not cal.is_business_day(Date(20, 7, 2026))

    # Mountain Day
    assert not cal.is_business_day(Date(11, 8, 2026))

    # Respect for the Aged Day
    assert not cal.is_business_day(Date(21, 9, 2026))

    # Citizens' holiday
    assert not cal.is_business_day(Date(22, 9, 2026))

    # Autumnal Equinox
    assert not cal.is_business_day(Date(23, 9, 2026))

    # Sports Day
    assert not cal.is_business_day(Date(12, 10, 2026))

    # Dec 31 is not a statutory Japanese public holiday
    assert cal.is_business_day(Date(31, 12, 2026))

def test_japan_emperors_birthday_history():

    cal = Calendar(CalendarTypes.JAPAN)

    # Previous Emperor
    assert not cal.is_business_day(Date(23, 12, 2018))

    # No Emperor's Birthday in 2019
    assert cal.is_business_day(Date(23, 12, 2019))

    # Current Emperor
    assert not cal.is_business_day(Date(23, 2, 2026))

def test_us_government_securities():

    cal = Calendar(
        CalendarTypes.US_GOVERNMENT_SECURITIES
    )
    
    # MLK
    assert not cal.is_business_day(
        Date(19, 1, 2026)
    )
    
    # Presidents Day
    assert not cal.is_business_day(
        Date(16, 2, 2026)
    )
    
    # Good Friday / no SOFR in 2026
    assert not cal.is_business_day(
        Date(3, 4, 2026)
    )
    
    # Memorial Day
    assert not cal.is_business_day(
        Date(25, 5, 2026)
    )
    
    # Juneteenth
    assert not cal.is_business_day(
        Date(19, 6, 2026)
    )
    
    # July 3 special 2026 closure / no SOFR
    assert not cal.is_business_day(
        Date(3, 7, 2026)
    )
    
    # July 4 is Saturday anyway
    assert not cal.is_business_day(
        Date(4, 7, 2026)
    )
    
    # Labor Day
    assert not cal.is_business_day(
        Date(7, 9, 2026)
    )
    
    # Columbus Day
    assert not cal.is_business_day(
        Date(12, 10, 2026)
    )
    
    # Veterans Day
    assert not cal.is_business_day(
        Date(11, 11, 2026)
    )
    
    # Thanksgiving
    assert not cal.is_business_day(
        Date(26, 11, 2026)
    )
    
    # Christmas
    assert not cal.is_business_day(
        Date(25, 12, 2026)
    )
    
    # Ordinary Treasury / SOFR business day
    assert cal.is_business_day(
        Date(6, 7, 2026)
    )


def test_us_federal_reserve():

    fed = Calendar(CalendarTypes.US_FEDERAL_RESERVE)
    ust = Calendar(CalendarTypes.US_GOVERNMENT_SECURITIES)
    
    # Federal Reserve Banks / EFFR calendar:
    assert fed.is_business_day(Date(3, 7, 2026)) is True
    
    # Government securities / SOFR calendar:
    assert ust.is_business_day(Date(3, 7, 2026)) is False

    cal = Calendar(CalendarTypes.US_FEDERAL_RESERVE)
    
    # July 4, 2026 is Saturday.
    # Fed Banks are open Friday July 3.
    assert cal.is_business_day(Date(3, 7, 2026)) is True
    
    # July 4 itself is weekend anyway.
    assert cal.is_business_day(Date(4, 7, 2026)) is False
    
    # July 4, 2027 is Sunday -> Monday July 5 observed.
    assert cal.is_business_day(Date(5, 7, 2027)) is False

def test_sydney():

    cal = Calendar(CalendarTypes.SYDNEY)
    
    # Australia Day
    assert not cal.is_business_day(Date(26, 1, 2026))
    
    # Good Friday
    assert not cal.is_business_day(Date(3, 4, 2026))
    
    # NSW additional ANZAC holiday
    assert not cal.is_business_day(Date(27, 4, 2026))
    
    # King's Birthday
    assert not cal.is_business_day(Date(8, 6, 2026))
    
    # NSW Bank Holiday
    assert not cal.is_business_day(Date(3, 8, 2026))
    
    # Labour Day
    assert not cal.is_business_day(Date(5, 10, 2026))
    
    # Christmas
    assert not cal.is_business_day(Date(25, 12, 2026))
    
    # Boxing Day additional holiday
    assert not cal.is_business_day(Date(28, 12, 2026))

def test_rits():

    rits = Calendar(CalendarTypes.AUSTRALIA_RITS)
    sydney = Calendar(CalendarTypes.SYDNEY)
    
    dt = Date(27, 4, 2026)
    
    assert rits.is_business_day(dt) is True
    assert sydney.is_business_day(dt) is False

    # NSW Bank Holiday
    assert Calendar(CalendarTypes.AUSTRALIA_RITS).is_business_day(
        Date(3, 8, 2026)
    ) is True
    
    assert Calendar(CalendarTypes.SYDNEY).is_business_day(
        Date(3, 8, 2026)
    ) is False
    
    # NSW Labour Day
    assert Calendar(CalendarTypes.AUSTRALIA_RITS).is_business_day(
        Date(5, 10, 2026)
    ) is True
    
    assert Calendar(CalendarTypes.SYDNEY).is_business_day(
        Date(5, 10, 2026)
    ) is False


def test_toronto():
    cal = Calendar(CalendarTypes.TORONTO)
    
    # Family Day
    assert not cal.is_business_day(Date(16, 2, 2026))
    
    # Good Friday
    assert not cal.is_business_day(Date(3, 4, 2026))
    
    # Victoria Day
    assert not cal.is_business_day(Date(18, 5, 2026))
    
    # Canada Day
    assert not cal.is_business_day(Date(1, 7, 2026))
    
    # Civic Holiday
    assert not cal.is_business_day(Date(3, 8, 2026))
    
    # Labour Day
    assert not cal.is_business_day(Date(7, 9, 2026))
    
    # Thanksgiving
    assert not cal.is_business_day(Date(12, 10, 2026))
    
    # Christmas
    assert not cal.is_business_day(Date(25, 12, 2026))

def test_zurich():

    cal = Calendar(CalendarTypes.ZURICH)

    assert not cal.is_business_day(Date(1, 1, 2026))
    assert not cal.is_business_day(Date(2, 1, 2026))
    assert not cal.is_business_day(Date(3, 4, 2026))
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(1, 5, 2026))
    assert not cal.is_business_day(Date(14, 5, 2026))
    assert not cal.is_business_day(Date(25, 5, 2026))
    assert not cal.is_business_day(Date(1, 8, 2026))
    assert not cal.is_business_day(Date(25, 12, 2026))
    assert not cal.is_business_day(Date(26, 12, 2026))

def test_hong_kong():

    cal = Calendar(CalendarTypes.HONG_KONG)
    
    assert not cal.is_business_day(Date(1, 1, 2026))
    
    # Lunar New Year
    assert not cal.is_business_day(Date(17, 2, 2026))
    assert not cal.is_business_day(Date(18, 2, 2026))
    assert not cal.is_business_day(Date(19, 2, 2026))
    
    # Good Friday
    assert not cal.is_business_day(Date(3, 4, 2026))
    
    # Ching Ming / Easter substitution sequence
    assert not cal.is_business_day(Date(6, 4, 2026))
    assert not cal.is_business_day(Date(7, 4, 2026))
    
    # Buddha's Birthday substitute
    assert not cal.is_business_day(Date(25, 5, 2026))
    
    # Dragon Boat / Tuen Ng
    assert not cal.is_business_day(Date(19, 6, 2026))
    
    # HKSAR Establishment Day
    assert not cal.is_business_day(Date(1, 7, 2026))
    
    # National Day
    assert not cal.is_business_day(Date(1, 10, 2026))
    
    # Chung Yeung substitute
    assert not cal.is_business_day(Date(19, 10, 2026))
    
    # Christmas
    assert not cal.is_business_day(Date(25, 12, 2026))

def test_tokyo():

    cal = Calendar(CalendarTypes.TOKYO)

    # BOJ New Year closures
    assert not cal.is_business_day(Date(2, 1, 2026))
    assert cal.is_holiday(Date(2, 1, 2026))

    # Jan 3 is also a BOJ holiday, although in 2026 it is Saturday.
    assert not cal.is_business_day(Date(3, 1, 2026))
    assert cal.is_holiday(Date(3, 1, 2026))

    # Japanese statutory holidays are also Tokyo holidays
    assert not cal.is_business_day(Date(20, 3, 2026))
    assert not cal.is_business_day(Date(6, 5, 2026))
    assert not cal.is_business_day(Date(20, 7, 2026))
    assert not cal.is_business_day(Date(11, 8, 2026))
    assert not cal.is_business_day(Date(21, 9, 2026))
    assert not cal.is_business_day(Date(22, 9, 2026))
    assert not cal.is_business_day(Date(23, 9, 2026))
    assert not cal.is_business_day(Date(12, 10, 2026))

    # BOJ year-end closure
    assert not cal.is_business_day(Date(31, 12, 2026))
    assert cal.is_holiday(Date(31, 12, 2026))


def test_aud_fx_settlement():
    fx = Calendar(CalendarTypes.AUD_FX_SETTLEMENT)
    syd = Calendar(CalendarTypes.SYDNEY)
    rits = Calendar(CalendarTypes.AUSTRALIA_RITS)

    # NSW-only holiday in 2026:
    dt = Date(27, 4, 2026)

    assert not fx.is_business_day(dt)
    assert not syd.is_business_day(dt)
    assert rits.is_business_day(dt)

def test_new_zealand_matariki_range():

    cal = Calendar(CalendarTypes.NEW_ZEALAND)

    # Last currently legislated Matariki date
    assert not cal.is_business_day(Date(21, 6, 2052))

    # Calendar should remain operational after the current
    # published Matariki schedule ends.
    dt = Date(15, 7, 2053)

    assert isinstance(cal.is_business_day(dt), bool)


def test_singapore():
    cal = Calendar(CalendarTypes.SINGAPORE)
    
    assert not cal.is_business_day(Date(1, 1, 2026))
    
    # Chinese New Year
    assert not cal.is_business_day(Date(17, 2, 2026))
    assert not cal.is_business_day(Date(18, 2, 2026))
    
    # Good Friday
    assert not cal.is_business_day(Date(3, 4, 2026))
    
    # Labour Day
    assert not cal.is_business_day(Date(1, 5, 2026))
    
    # Hari Raya Haji
    assert not cal.is_business_day(Date(27, 5, 2026))
    
    # Vesak observed Monday
    assert not cal.is_business_day(Date(1, 6, 2026))
    
    # National Day observed Monday
    assert not cal.is_business_day(Date(10, 8, 2026))
    
    # Deepavali observed Monday
    assert not cal.is_business_day(Date(9, 11, 2026))
    
    # Christmas
    assert not cal.is_business_day(Date(25, 12, 2026))

def test_new_york_historical_juneteenth():

    cal = Calendar(CalendarTypes.NEW_YORK)

    # Before Juneteenth became a federal holiday
    assert cal.is_business_day(Date(19, 6, 2020))

    # After introduction
    assert not cal.is_business_day(Date(19, 6, 2023))

###############################################################################

###############################################################################

test_calendar()

# Market calendar tests
test_weekend()
test_australia()
test_canada()
test_france()
test_germany()
test_italy()
test_new_zealand()
test_norway()
test_sweden()
test_switzerland()

test_target()
test_target_not_national_holidays()

test_united_states()
test_new_york()
test_us_government_securities()
test_us_federal_reserve()

test_london()
test_united_kingdom()

test_japan()
test_tokyo()

test_sydney()
test_rits()
test_aud_fx_settlement()

test_toronto()
test_zurich()
test_hong_kong()
test_singapore()

test_japan_emperors_birthday_history()
test_new_zealand_matariki_range()
test_new_york_historical_juneteenth()

# Calendar engine / joint-calendar tests
test_fast_adjust()
test_joint_calendar_holiday_is_union()
test_joint_business_day_is_intersection_over_full_years()
test_joint_adjust_skips_holidays_of_both_markets()
test_joint_fast_adjust_matches_adjust()
test_joint_add_business_days_round_trip()
test_joint_holiday_list_is_union_of_lists()
test_joint_single_entry_and_duplicates_collapse()
test_joint_str_and_repr()

test_cases.compare_test_cases()