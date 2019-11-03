# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.market.curves.FinFlatCurve import FinFlatCurve, FinCompoundingMethods
from financepy.finutils.FinDate import FinDate
import numpy as np
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinFlatCurve():

    curveDate = FinDate(2019, 2, 2)
    times = np.linspace(0.0, 1.0, 5)

    testCases.header("COMPOUNDING", "DFS")
    compounding = FinCompoundingMethods.CONTINUOUS
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)

    compounding = FinCompoundingMethods.ANNUAL
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)

    compounding = FinCompoundingMethods.SEMI_ANNUAL
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)

    compounding = FinCompoundingMethods.QUARTERLY
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)

    compounding = FinCompoundingMethods.MONTHLY
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)

    compounding = FinCompoundingMethods.MONEY_MARKET
    flatCurve = FinFlatCurve(curveDate, 0.05, compounding)
    dfs = flatCurve.df(times)
    testCases.print(compounding, dfs)


test_FinFlatCurve()
testCases.compareTestCases()
