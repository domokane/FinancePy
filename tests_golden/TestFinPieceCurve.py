###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")

# import numpy as np


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################
#
# def test_FinPieceCurve():
#
#    times = np.linspace(0.0, 1.0, 5)
#    values = np.ones(5) * 0.05
#
#    flat_curve = FinPieceCurve(times,values)
#
#    dfs = flat_curve.df(times, 0)
#    testCases.print(dfs)
#    dfs = flat_curve.df(times, 1)
#    testCases.print(dfs)
#    dfs = flat_curve.df(times, 2)
#    testCases.print(dfs)
#    dfs = flat_curve.df(times, -1)
#    testCases.print(dfs)

##########################################################################


# test_FinPieceCurve()
testCases.compareTestCases()
