###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinDayCount():

    testCases.header("DAY_COUNT_METHOD", "START", "END", "ALPHA")

    finFreq = FrequencyTypes.ANNUAL

    for day_count_method in DayCountTypes:

        start_date = Date(1, 1, 2019)
        next_date = start_date
        numDays = 20
        day_count = DayCount(day_count_method)

        for _ in range(0, numDays):
            next_date = next_date.add_days(7)
            dcf = day_count.year_frac(
                start_date, next_date, next_date, finFreq)

            testCases.print(
                str(day_count_method),
                str(start_date),
                str(next_date),
                dcf[0])


test_FinDayCount()
testCases.compareTestCases()
