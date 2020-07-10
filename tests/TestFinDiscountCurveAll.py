# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.market.curves.FinInterpolate import FinInterpMethods

from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros

import numpy as np
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
# TODO: Add other discount curves
###############################################################################


def test_FinDiscountCurves():

    # Create a curve from times and discount factors
    startDate = FinDate(1, 1, 2018)
    years = np.linspace(0, 10, 11)
    rate = 0.05
    values = np.exp(-rate * years)
    dates = startDate.addYears(years)

    curve = FinDiscountCurve(startDate,
                             dates,
                             values,
                             FinInterpMethods.FLAT_FORWARDS)

    testCases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    date1 = FinDate(1, 6, 2018)
    date2 = FinDate(1, 1, 2019)
    basisType = FinDayCountTypes.ACT_365_ISDA

    years = np.linspace(0, 10, 11)
    dates = startDate.addYears(years)

    for dt in dates:
        df = curve.df(dt)
        zeroRate = curve.zeroRate(dt, FinFrequencyTypes.ANNUAL)
        fwd = curve.fwd(dt)
        fwdRate = curve.fwdRate(date1, date2, basisType)
        q = curve.survProb(dt)
        testCases.print(dt, df, zeroRate, fwd, fwdRate, q)

    ###########################################################################
    # Curve built from a single rate so is flat
    ###########################################################################
    flatRate = 0.05
    curve = FinDiscountCurveFlat(startDate,
                                 flatRate)

    testCases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    date1 = FinDate(1, 6, 2018)
    date2 = FinDate(1, 1, 2019)
    basisType = FinDayCountTypes.ACT_365_ISDA

    for dt in dates:
        df = curve.df(dt)
        zeroRate = curve.zeroRate(dt, FinFrequencyTypes.ANNUAL)
        fwd = curve.fwd(dt)
        fwdRate = curve.fwdRate(date1, date2, basisType)
        q = curve.survProb(dt)
        testCases.print(dt, df, zeroRate, fwd, fwdRate, q)

    ###########################################################################
    # Curve Built from Zero Rates
    ###########################################################################

    flatRate = 0.05
    times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    dates = startDate.addYears(times)
    zeroRates = []
    for i in range(0, len(times)):
        r = 0.05 + 0.0020*i
        zeroRates.append(r)

    curve = FinDiscountCurveZeros(startDate,
                                  dates,
                                  zeroRates)

    testCases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    date1 = FinDate(1, 6, 2018)
    date2 = FinDate(1, 1, 2019)
    basisType = FinDayCountTypes.ACT_365_ISDA

    for dt in dates:
        df = curve.df(dt)
        zeroRate = curve.zeroRate(dt, FinFrequencyTypes.ANNUAL)
        fwd = curve.fwd(dt)
        fwdRate = curve.fwdRate(date1, date2, basisType)
        q = curve.survProb(dt)
        testCases.print(dt, df, zeroRate, fwd, fwdRate, q)

###############################################################################


test_FinDiscountCurves()
testCases.compareTestCases()
