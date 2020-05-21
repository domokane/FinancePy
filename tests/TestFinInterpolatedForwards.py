# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinInterpolate import interpolate, FinInterpMethods
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinFlatCurve import FinFlatCurve

import numpy as np
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = True

def test_FinInterpolatedForwards():

    import matplotlib.pyplot as plt

    tValues = np.array([0.0, 3.0, 5.0, 10.0])
    rValues = np.array([0.04, 0.07, 0.08, 0.09])
    dfValues = np.exp(-tValues*rValues)
    tInterpValues = np.linspace(0.0, 12.0, 200)

    print(tValues)
    print(rValues)
    print(dfValues)

    curveDate = FinDate(2019,1,3)
    for method in FinInterpMethods:

        discountCurve = FinDiscountCurve(curveDate, tValues, dfValues, method)
        dfInterpValues = discountCurve.df(tInterpValues)
        fwdInterpValues = discountCurve.fwd(tInterpValues)
        zeroInterpValues = discountCurve.zeroRate(tInterpValues)

        if PLOT_GRAPHS:
            plt.figure(figsize=(8, 6))
            plt.plot(tValues, dfValues, 'o', color='g', label="DFS:")
            plt.plot(tInterpValues, dfInterpValues, color='r', label="DF:" + str(method))
            plt.legend()
            plt.figure(figsize=(8, 6))
            plt.plot(tInterpValues, fwdInterpValues, color='r', label="FWD:" + str(method))
            plt.plot(tInterpValues, zeroInterpValues, color='b', label="ZERO:" + str(method))
            plt.plot(tValues, rValues, 'o', color='g',  label="ZERO RATES")
            plt.legend()

test_FinInterpolatedForwards()
testCases.compareTestCases()
