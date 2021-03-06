###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinInterpolator import FinInterpTypes
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinInterpolatedForwards():

    import matplotlib.pyplot as plt

    tValues = np.array([0.0, 3.0, 5.0, 10.0])
    rValues = np.array([0.04, 0.07, 0.08, 0.09])
    dfValues = np.exp(-tValues*rValues)
    tInterpValues = np.linspace(0.0, 12.0, 49)

    curveDate = FinDate(1, 1, 2019)

    tDates = curveDate.addYears(tValues)
    tInterpDates = curveDate.addYears(tInterpValues)

    for interpType in FinInterpTypes:

        discountCurve = FinDiscountCurve(curveDate, tDates, dfValues, interpType)
        dfInterpValues = discountCurve.df(tInterpDates)
        fwdInterpValues = discountCurve.fwd(tInterpDates)
        zeroInterpValues = discountCurve.zeroRate(tInterpDates)

        if PLOT_GRAPHS:
            plt.figure(figsize=(8, 6))
            plt.plot(tValues, dfValues, 'o', color='g', label="DFS:")
            plt.plot(tInterpValues, dfInterpValues, color='r',
                     label="DF:" + str(interpType))
            plt.legend()
            plt.figure(figsize=(8, 6))
            plt.plot(tInterpValues, fwdInterpValues, color='r',
                     label="FWD:" + str(interpType))
            plt.plot(tInterpValues, zeroInterpValues, color='b',
                     label="ZERO:" + str(interpType))
            plt.plot(tValues, rValues, 'o', color='g',  label="ZERO RATES")
            plt.legend()

###############################################################################


test_FinInterpolatedForwards()
testCases.compareTestCases()
