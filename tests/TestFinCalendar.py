###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDate import setDateFormatType, FinDateFormatTypes
from financepy.finutils.FinCalendar import FinCalendar, FinCalendarTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinCalendar():

    setDateFormatType(FinDateFormatTypes.US_LONGEST)
    endDate = FinDate(31, 12, 2030)

    for calendarType in FinCalendarTypes:

        testCases.banner("================================")
        testCases.banner("================================")

        testCases.header("CALENDAR", "HOLIDAY")
        testCases.print("STARTING", calendarType)

        cal = FinCalendar(calendarType)
        nextDate = FinDate(31, 12, 2020)

        while nextDate < endDate:
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


test_FinCalendar()
testCases.compareTestCases()
