##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_Amount():

    testCases.header("LABEL", "AMOUNT")
    x = Amount(101000.232)

    testCases.print("Amount", x)

    x = Amount(101000.232, CurrencyTypes.CAD)

    testCases.print("Amount", x)

###############################################################################


test_Amount()

testCases.compareTestCases()

