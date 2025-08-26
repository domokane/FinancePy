# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.calendar import Calendar, CalendarTypes
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date


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


########################################################################################

test_calendar()
test_cases.compare_test_cases()
