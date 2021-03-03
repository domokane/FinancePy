##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinAmount import FinAmount
from financepy.finutils.FinCurrency import FinCurrencyTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################

def test_FinAmount():

    testCases.header("LABEL", "AMOUNT")
    x = FinAmount(101000.232)

    testCases.print("Amount", x)

    x = FinAmount(101000.232, FinCurrencyTypes.CAD)

    testCases.print("Amount", x)

###############################################################################

test_FinAmount()

testCases.compareTestCases()
