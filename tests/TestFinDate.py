# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinDate():

    startDate = FinDate(2018, 1, 1)

    testCases.header("DATE", "MONTHS", "CDS DATE")

    for numMonths in range(0, 120):
        nextCDSDate = startDate.nextCDSDate(numMonths)
        testCases.print(str(startDate), numMonths, str(nextCDSDate))

    startDate = FinDate(2018, 1, 1)

    testCases.header("STARTDATE", "MONTHS", "CDS DATE")

    for numMonths in range(0, 365):
        startDate = startDate.addDays(1)
        nextIMMDate = startDate.nextIMMDate()
        testCases.print(numMonths, str(startDate), str(nextIMMDate))


def test_FinDateTenors():

    startDate = FinDate(2018, 2, 23)

    tenor = "5d"
    print(tenor, startDate.addTenor(tenor))

    tenor = "7D"
    print(tenor, startDate.addTenor(tenor))

    tenor = "1W"
    print(tenor, startDate.addTenor(tenor))

    tenor = "4W"
    print(tenor, startDate.addTenor(tenor))

    tenor = "1M"
    print(tenor, startDate.addTenor(tenor))

    tenor = "24M"
    print(tenor, startDate.addTenor(tenor))

    tenor = "2Y"
    print(tenor, startDate.addTenor(tenor))

    tenor = "10y"
    print(tenor, startDate.addTenor(tenor))

    tenor = "0m"
    print(tenor, startDate.addTenor(tenor))

    tenor = "20Y"
    print(tenor, startDate.addTenor(tenor))


test_FinDate()
test_FinDateTenors()

testCases.compareTestCases()
