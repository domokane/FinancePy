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
from financepy.products.bonds.FinMortgage import FinMortgage, FinMortgageType
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

def test_FinMortgage():

    principal = 130000
    startDate = FinDate(23, 2, 2018)
    endDate = startDate.addTenor("10Y")
    mortgage = FinMortgage(startDate, endDate, principal)

    rate = 0.035
    mortgage.generateFlows(rate, FinMortgageType.REPAYMENT)
    mortgage.print()

    mortgage.generateFlows(rate, FinMortgageType.INTEREST_ONLY)
    mortgage.print()

###############################################################################

test_FinMortgage()
testCases.compareTestCases()
