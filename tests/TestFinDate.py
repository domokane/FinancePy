# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinDate import FinDate, dateRange
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


def test_FinDateHelpers():
    startDate = FinDate(2010, 1, 1)
    
    testCases.header("Tenor", "Dates")
    
    endDate = startDate.addDays(3)
    tenor = "Default"
    testCases.print(tenor, dateRange(startDate, endDate))
    
    endDate = startDate.addDays(20)
    tenor = "1W"
    testCases.print(tenor, dateRange(startDate, endDate, tenor))
    
    tenor = "7D"
    testCases.print(tenor, dateRange(startDate, endDate, tenor))
    
    testCases.header("Case", "Dates")
    
    case = "Same startDate"
    testCases.print(case, dateRange(startDate, startDate))
    case = "startDate before endDate"
    testCases.print(case, dateRange(endDate, startDate))


test_FinDate()
test_FinDateTenors()
test_FinDateHelpers()

testCases.compareTestCases()
