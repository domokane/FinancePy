# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from FinTestCases import FinTestCases, global_test_case_mode
import sys

sys.path.append("..")

# import numpy as np


test_cases = FinTestCases(__file__, global_test_case_mode)

# def test_FinPieceCurve():
#    times = np.linspace(0.0, 1.0, 5)
#    values = np.ones(5) * 0.05
#    flat_curve = FinPieceCurve(times,values)
#    dfs = flat_curve.df(times, 0)
#    test_cases.print(dfs)
#    dfs = flat_curve.df(times, 1)
#    test_cases.print(dfs)
#    dfs = flat_curve.df(times, 2)
#    test_cases.print(dfs)
#    dfs = flat_curve.df(times, -1)
#    test_cases.print(dfs)



# test_FinPieceCurve()
test_cases.compare_test_cases()
