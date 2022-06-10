###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.date import Date


freq_type = FrequencyTypes.SEMI_ANNUAL
date_gen_rule_type = DateGenRuleTypes.BACKWARD


def test_date_adjust_no_adj():
    start_date = Date(28, 2, 2008)
    end_date = Date(28, 2, 2011)

    calendar_type = CalendarTypes.NONE
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    assert schedule._adjusted_dates == [
        Date(28, 2, 2008), Date(28, 8, 2008), Date(
            28, 2, 2009), Date(28, 8, 2009),
        Date(28, 2, 2010), Date(28, 8, 2010), Date(28, 2, 2011)]


def test_date_adjust_noweekend_following():
    start_date = Date(28, 2, 2008)
    end_date = Date(28, 2, 2011)

    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    assert schedule._adjusted_dates == [
        Date(28, 2, 2008), Date(28, 8, 2008), Date(
            2, 3, 2009), Date(28, 8, 2009),
        Date(1, 3, 2010), Date(30, 8, 2010), Date(28, 2, 2011)]


def test_date_adjust_noweekend_modfollowing():
    start_date = Date(28, 2, 2008)
    end_date = Date(28, 2, 2011)

    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    assert schedule._adjusted_dates == [
        Date(28, 2, 2008), Date(28, 8, 2008), Date(
            27, 2, 2009), Date(28, 8, 2009),
        Date(26, 2, 2010), Date(30, 8, 2010), Date(28, 2, 2011)]


def test_date_adjust_noweekend_usholidays_modfollowing():
    start_date = Date(4, 7, 2008)
    end_date = Date(4, 7, 2011)

    calendar_type = CalendarTypes.UNITED_STATES
    bus_day_adjust_type = BusDayAdjustTypes.MODIFIED_FOLLOWING

    schedule = Schedule(start_date,
                        end_date,
                        freq_type,
                        calendar_type,
                        bus_day_adjust_type,
                        date_gen_rule_type)

    assert schedule._adjusted_dates == [
        Date(4, 7, 2008), Date(5, 1, 2009), Date(6, 7, 2009), Date(4, 1, 2010),
        Date(6, 7, 2010), Date(4, 1, 2011), Date(5, 7, 2011)]
