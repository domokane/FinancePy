###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date

test_cases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinDayCount():

    test_cases.header("DAY_COUNT_METHOD", "START", "END", "ALPHA")

    finFreq = FrequencyTypes.ANNUAL

    for day_count_method in DayCountTypes:

        start_dt = Date(1, 1, 2019)
        next_dt = start_dt
        num_days = 20
        day_count = DayCount(day_count_method)

        for _ in range(0, num_days):
            next_dt = next_dt.add_days(7)
            dcf = day_count.year_frac(
                start_dt, next_dt, next_dt, finFreq)

            test_cases.print(
                str(day_count_method),
                str(start_dt),
                str(next_dt),
                dcf[0])


test_FinDayCount()
test_cases.compareTestCases()
