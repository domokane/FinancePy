###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.date import setDateFormatType, FinDateFormatTypes
from financepy.utils.calendar import Calendar, CalendarTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_Calendar():

    setDateFormatType(FinDateFormatTypes.US_LONGEST)
    end_date = Date(31, 12, 2030)

    for calendar_type in CalendarTypes:

        testCases.banner("================================")
        testCases.banner("================================")

        testCases.header("CALENDAR", "HOLIDAY")
        testCases.print("STARTING", calendar_type)

        cal = Calendar(calendar_type)
        next_date = Date(31, 12, 2020)

        while next_date < end_date:
            next_date = next_date.addDays(1)
            
            if next_date._d == 1 and next_date._m == 1:
                testCases.banner("================================")
#                print("=========================")

            isHolidayDay = cal.isHoliday(next_date)
            if isHolidayDay is True:
                testCases.print(cal, next_date)
#                print(cal, next_date)

    setDateFormatType(FinDateFormatTypes.US_LONG)

###############################################################################


test_Calendar()
testCases.compareTestCases()
