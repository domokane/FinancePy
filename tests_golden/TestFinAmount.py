##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount
from FinTestCases import FinTestCases, globalTestCaseMode

import sys
sys.path.append("..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_amount():

    testCases.header("LABEL", "AMOUNT")
    x = Amount(101000.232)

    testCases.print("Amount", x)

    x = Amount(101000.232, CurrencyTypes.CAD)

    testCases.print("Amount", x)

###############################################################################


test_amount()

testCases.compareTestCases()
