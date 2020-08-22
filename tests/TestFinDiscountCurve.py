###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.market.curves.FinInterpolate import FinInterpTypes

from financepy.market.curves.FinDiscountCurve import FinDiscountCurve

from financepy.finutils.FinMath import scale

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
# TODO: Add other discount curves
###############################################################################

PLOT_GRAPHS = False


def test_FinDiscountCurve():

    # Create a curve from times and discount factors
    start = FinDate(1, 1, 2018)
    years = np.linspace(0, 10, 6)
    rate = 0.05 + 0.005*years - 0.0003*years*years
    dfs = np.exp(-rate * years)
    dates = start.addYears(years)

    curve = FinDiscountCurve(start, dates, dfs, FinInterpTypes.FLAT_FORWARDS)

    testCases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    plotYears = np.linspace(0, 12, 12*12+1)
    plotDates = start.addYears(plotYears)

    # Examine dependency of curve on compounding rate
    zeroRates_A = curve.zeroRate(plotDates, FinFrequencyTypes.ANNUAL)
    zeroRates_S = curve.zeroRate(plotDates, FinFrequencyTypes.SEMI_ANNUAL)
    zeroRates_Q = curve.zeroRate(plotDates, FinFrequencyTypes.QUARTERLY)
    zeroRates_M = curve.zeroRate(plotDates, FinFrequencyTypes.MONTHLY)
    zeroRates_C = curve.zeroRate(plotDates, FinFrequencyTypes.CONTINUOUS)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(plotYears, scale(zeroRates_A, 100), label='A')
        plt.plot(plotYears, scale(zeroRates_S, 100), label='S')
        plt.plot(plotYears, scale(zeroRates_Q, 100), label='Q')
        plt.plot(plotYears, scale(zeroRates_M, 100), label='M')
        plt.plot(plotYears, scale(zeroRates_C, 100), label='C')
        plt.ylim((5, 8))

        plt.title('Discount Curves')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    # Examine dependency of fwd curve on the interpolation scheme

    for interp in FinInterpTypes:

        curve = FinDiscountCurve(start, dates, dfs, interp)
        fwdRates = curve.fwd(plotDates)
        zeroRates = curve.zeroRate(plotDates, FinFrequencyTypes.ANNUAL)
        parRates = curve.parRate(plotDates, FinFrequencyTypes.ANNUAL)

        if PLOT_GRAPHS:
            plt.figure(figsize=(6, 4))
            plt.plot(plotYears, scale(fwdRates, 100), label='FWD RATES')
            plt.plot(plotYears, scale(zeroRates, 100), label='ZERO RATES')
            plt.plot(plotYears, scale(parRates, 100), label='PAR RATES')
            plt.ylim((3.0, 8.5))

            plt.title('Forward Curves using ' + str(interp))
            plt.xlabel('Time (years)')
            plt.ylabel('Fwd Rate (%)')
            plt.legend(loc='lower right', frameon=False)

###############################################################################


test_FinDiscountCurve()
testCases.compareTestCases()
