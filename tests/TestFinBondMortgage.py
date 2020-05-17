# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:12 2019

@author: Dominic
"""

import sys

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinCalendar import FinCalendarTypes
from financepy.products.bonds.FinBondMortgage import FinBondMortgage, FinBondMortgageType
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinBondMortgage():

    principal = 130000
    startDate = FinDate(23, 2, 2018)
    endDate = startDate.addTenor("10Y")
    mortgage = FinBondMortgage(startDate, endDate, principal)

    rate = 0.035
    mortgage.generateFlows(rate, FinBondMortgageType.REPAYMENT)

    numFlows = len(mortgage._schedule._adjustedDates)

    testCases.header("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING",
                     "TOTAL")

    for i in range(0, numFlows):
        testCases.print(mortgage._schedule._adjustedDates[i],
                        mortgage._interestFlows[i],
                        mortgage._principalFlows[i],
                        mortgage._principalRemaining[i],
                        mortgage._totalFlows[i])

    mortgage.generateFlows(rate, FinBondMortgageType.INTEREST_ONLY)

    testCases.header("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING",
                     "TOTAL")

    for i in range(0, numFlows):
        testCases.print(mortgage._schedule._adjustedDates[i],
                        mortgage._interestFlows[i],
                        mortgage._principalFlows[i],
                        mortgage._principalRemaining[i],
                        mortgage._totalFlows[i])
###############################################################################


test_FinBondMortgage()
testCases.compareTestCases()
