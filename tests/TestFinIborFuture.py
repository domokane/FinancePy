###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.rates.FinIborFuture import FinIborFuture
from financepy.finutils.FinDate import *

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

setDateFormatType(FinDateFormatTypes.UK_LONG)

###############################################################################


def test_FinIborFuture():

    todayDate = FinDate(5, 5, 2020)

    testCases.header("VALUES")

    for i in range(1, 12):
        fut = FinIborFuture(todayDate, i, "3M")
        testCases.print(fut)

        fra = fut.toFRA(0.020, 0.0)
        testCases.print(fra)

###############################################################################


test_FinIborFuture()
testCases.compareTestCases()
