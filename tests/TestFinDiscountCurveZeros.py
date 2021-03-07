###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinInterpolator import FinInterpTypes
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinDiscountCurveZeros():

    startDate = FinDate(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    dates = startDate.addYears(times)
    zeroRates = np.linspace(5.0, 6.0, 10)/100
    freqType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA

    curve = FinDiscountCurveZeros(startDate,
                                  dates,
                                  zeroRates,
                                  freqType,
                                  dayCountType,
                                  FinInterpTypes.FLAT_FWD_RATES)

    testCases.header("T", "DF")

    years = np.linspace(0, 10, 21)
    dates = startDate.addYears(years)
    for dt in dates:
        df = curve.df(dt)
        testCases.print(dt, df)

#    print(curve)

###############################################################################

    numRepeats = 100

    start = time.time()

    for i in range(0, numRepeats):
        freqType = FinFrequencyTypes.ANNUAL
        dayCountType = FinDayCountTypes.ACT_ACT_ISDA

        dates = [FinDate(14, 6, 2016), FinDate(14, 9, 2016),
                 FinDate(14, 12, 2016), FinDate(14, 6, 2017),
                 FinDate(14, 6, 2019), FinDate(14, 6, 2021),
                 FinDate(15, 6, 2026), FinDate(16, 6, 2031),
                 FinDate(16, 6, 2036), FinDate(14, 6, 2046)]

        zeroRates = [0.000000, 0.006616, 0.007049, 0.007795,
                     0.009599, 0.011203, 0.015068, 0.017583,
                     0.018998, 0.020080]

        startDate = dates[0]

        curve = FinDiscountCurveZeros(startDate,
                                      dates,
                                      zeroRates,
                                      freqType,
                                      dayCountType,
                                      FinInterpTypes.FLAT_FWD_RATES)

    end = time.time()
    period = end - start
#    print("Time taken:", period)

#    print(curve)

###############################################################################


test_FinDiscountCurveZeros()
testCases.compareTestCases()
