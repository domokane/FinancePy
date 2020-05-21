# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""
import time

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinInterpolate import FinInterpMethods
from financepy.market.curves.FinZeroCurve import FinZeroCurve
import numpy as np
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinDiscountCurve():

    startDate = FinDate(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    zeroRates = np.linspace(5.0, 6.0, 10)/100
    freqType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_ACT_ISDA

    curve = FinZeroCurve(startDate,
                         times,
                         zeroRates,
                         freqType,
                         dayCountType,
                         FinInterpMethods.FLAT_FORWARDS)

    testCases.header("T", "DF")

    for t in np.linspace(0, 10, 21):
        df = curve.df(t)
        testCases.print(t, df)

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

        curve = FinZeroCurve(startDate,
                             dates,
                             zeroRates,
                             freqType,
                             dayCountType,
                             FinInterpMethods.FLAT_FORWARDS)

#    end = time.time()
#    period = end - start
#    print("Time taken:", period)

#    print(curve)

###############################################################################


test_FinDiscountCurve()
testCases.compareTestCases()
