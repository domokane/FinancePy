# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np

from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.finutils.FinInterpolate import FinInterpMethods
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinDiscountCurve():

    startDate = FinDate(2018,1,1)
    times = np.linspace(0,10.0,10)
    rate= 0.05
    values = np.exp(-rate * times)
    
    curve = FinDiscountCurve(startDate,times,values, FinInterpMethods.FLAT_FORWARDS)

    testCases.header("T","DF")
    
    for t in np.linspace(0,10,21):
        df = curve.df(t)
        testCases.print(t,df)

test_FinDiscountCurve()
testCases.compareTestCases()
