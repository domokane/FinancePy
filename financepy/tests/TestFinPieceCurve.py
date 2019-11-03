# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..//..")

#import numpy as np

#from market.curves.FinPieceCurve import FinPieceCurve

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################

# def test_FinPieceCurve():
#
#    times = np.linspace(0.0,1.0,5)
#    values = np.ones(5) * 0.05
#
##    flatCurve = FinPieceCurve(times,values)
#
#    dfs = flatCurve.df(times,0)
#    testCases.print(dfs)
#    dfs = flatCurve.df(times,1)
#    testCases.print(dfs)
#    dfs = flatCurve.df(times,2)
#    testCases.print(dfs)
#    dfs = flatCurve.df(times,-1)
#    testCases.print(dfs)

##########################################################################

# NEED TO FIX CODE
# test_FinPieceCurve()
# testCases.compareTestCases()
