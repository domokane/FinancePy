###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate

from financepy.products.bonds.FinBondMortgage import FinBondMortgage
from financepy.products.bonds.FinBondMortgage import FinBondMortgageTypes
from financepy.products.rates.FinIborSingleCurve import FinIborSingleCurve

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)


###############################################################################


def test_FinBondMortgage():

    principal = 130000
    startDate = FinDate(23, 2, 2018)
    endDate = startDate.addTenor("10Y")
    mortgage = FinBondMortgage(startDate, endDate, principal)

    rate = 0.035
    mortgage.generateFlows(rate, FinBondMortgageTypes.REPAYMENT)

    numFlows = len(mortgage._schedule._adjustedDates)

    testCases.header("PAYMENT DATE", "INTEREST", "PRINCIPAL", "OUTSTANDING",
                     "TOTAL")

    for i in range(0, numFlows):
        testCases.print(mortgage._schedule._adjustedDates[i],
                        mortgage._interestFlows[i],
                        mortgage._principalFlows[i],
                        mortgage._principalRemaining[i],
                        mortgage._totalFlows[i])

    mortgage.generateFlows(rate, FinBondMortgageTypes.INTEREST_ONLY)

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
