###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.Date import Date
from financepy.utils.DayCount import DayCount, FinDayCountTypes
from financepy.utils.Frequency import FinFrequencyTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinDayCount():

    testCases.header("DAY_COUNT_METHOD", "START", "END", "ALPHA")

    finFreq = FinFrequencyTypes.ANNUAL

    for dayCountMethod in FinDayCountTypes:

        start_date = Date(1, 1, 2019)
        nextDate = start_date
        numDays = 20
        dayCount = DayCount(dayCountMethod)

        for _ in range(0, numDays):
            nextDate = nextDate.addDays(7)
            dcf = dayCount.year_frac(start_date, nextDate, nextDate, finFreq)

            testCases.print(
                str(dayCountMethod),
                str(start_date),
                str(nextDate),
                dcf[0])


test_FinDayCount()
testCases.compareTestCases()
