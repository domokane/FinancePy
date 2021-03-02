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
        nextDate = Date(31, 12, 2020)

        while nextDate < end_date:
            nextDate = nextDate.addDays(1)
            
            if nextDate._d == 1 and nextDate._m == 1:
                testCases.banner("================================")
#                print("=========================")

            isHolidayDay = cal.isHoliday(nextDate)
            if isHolidayDay is True:
                testCases.print(cal, nextDate)
#                print(cal, nextDate)

    setDateFormatType(FinDateFormatTypes.US_LONG)

###############################################################################


test_Calendar()
testCases.compareTestCases()
