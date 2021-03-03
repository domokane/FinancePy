###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.date import Date
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.frequency import FrequencyTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinDayCount():

    testCases.header("DAY_COUNT_METHOD", "START", "END", "ALPHA")

    finFreq = FrequencyTypes.ANNUAL

    for dayCountMethod in DayCountTypes:

        start_date = Date(1, 1, 2019)
        next_date = start_date
        numDays = 20
        dayCount = DayCount(dayCountMethod)

        for _ in range(0, numDays):
            next_date = next_date.addDays(7)
            dcf = dayCount.year_frac(start_date, next_date, next_date, finFreq)

            testCases.print(
                str(dayCountMethod),
                str(start_date),
                str(next_date),
                dcf[0])


test_FinDayCount()
testCases.compareTestCases()
