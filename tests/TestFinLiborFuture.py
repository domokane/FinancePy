# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.libor.FinLiborFuture import FinLiborFuture
from financepy.finutils.FinDate import FinDate

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinLiborFuture():

    todayDate = FinDate(5, 5, 2020)

    for i in range(1, 12):
        fut = FinLiborFuture(todayDate, i, "3M")
        print(fut)

        fra = fut.toFRA(0.020, 0.0)
        print(fra)

###############################################################################


test_FinLiborFuture()
testCases.compareTestCases()
