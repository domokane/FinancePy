###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np

import sys
sys.path.append("..")

from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.interpolator import FinInterpTypes
from financepy.market.curves.FinDiscountCurveZeros import DiscountCurveZeros

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinDiscountCurveZeros():

    start_date = Date(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    dates = start_date.addYears(times)
    zeroRates = np.linspace(5.0, 6.0, 10)/100
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.ACT_ACT_ISDA

    curve = DiscountCurveZeros(start_date,
                               dates,
                               zeroRates,
                               freq_type,
                               day_count_type,
                               FinInterpTypes.FLAT_FWD_RATES)

    testCases.header("T", "DF")

    years = np.linspace(0, 10, 21)
    dates = start_date.addYears(years)
    for dt in dates:
        df = curve.df(dt)
        testCases.print(dt, df)

#    print(curve)

###############################################################################

    numRepeats = 100

    start = time.time()

    for i in range(0, numRepeats):
        freq_type = FrequencyTypes.ANNUAL
        day_count_type = DayCountTypes.ACT_ACT_ISDA

        dates = [Date(14, 6, 2016), Date(14, 9, 2016),
                 Date(14, 12, 2016), Date(14, 6, 2017),
                 Date(14, 6, 2019), Date(14, 6, 2021),
                 Date(15, 6, 2026), Date(16, 6, 2031),
                 Date(16, 6, 2036), Date(14, 6, 2046)]

        zeroRates = [0.000000, 0.006616, 0.007049, 0.007795,
                     0.009599, 0.011203, 0.015068, 0.017583,
                     0.018998, 0.020080]

        start_date = dates[0]

        curve = DiscountCurveZeros(start_date,
                                   dates,
                                   zeroRates,
                                   freq_type,
                                   day_count_type,
                                   FinInterpTypes.FLAT_FWD_RATES)

    end = time.time()
    period = end - start
#    print("Time taken:", period)

#    print(curve)

###############################################################################


test_FinDiscountCurveZeros()
testCases.compareTestCases()
