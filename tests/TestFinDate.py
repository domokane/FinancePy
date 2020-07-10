# -*- coding: utf-8 -*-
"""
Created on Sat Feb 06 07:26:46 2016

@author: Dominic O'Kane
"""

import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate, dateRange
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


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

###############################################################################


def test_FinDateTenors():

    startDate = FinDate(2018, 2, 23)

    testCases.header("TENOR", "DATE")
    tenor = "5d"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "7D"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "1W"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "4W"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "1M"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "24M"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "2Y"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "10y"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "0m"
    testCases.print(tenor, startDate.addTenor(tenor))

    tenor = "20Y"
    testCases.print(tenor, startDate.addTenor(tenor))

###############################################################################


def test_FinDateRange():
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

###############################################################################


def test_FinDateAddMonths():
    startDate = FinDate(2010, 1, 1)

    testCases.header("Tenor", "Dates")

    months = [1, 3, 6, 9, 12, 24, 36, 48, 60]
    dates1 = startDate.addMonths(months)
    months = np.array(months)
    dates2 = startDate.addMonths(months)

    years = [1, 3, 5, 7, 10]
    dates3 = startDate.addYears(years)
    years = np.array(years)
    dates4 = startDate.addYears(years)
    print(dates4)

    years = [1.5, 3.25, 5.75, 7.25, 10.0]
    dates5 = startDate.addYears(years)
    years = np.array(years)
    dates6 = startDate.addYears(years)
    print(dates5)

###############################################################################


test_FinDate()
test_FinDateTenors()
test_FinDateRange()
test_FinDateAddMonths()

testCases.compareTestCases()
