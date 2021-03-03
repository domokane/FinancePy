###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityForward import FinEquityForward
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinLongShort
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityForward():

    valueDate = FinDate(13, 2, 2018)
    expiryDate = valueDate.addMonths(12)

    stockPrice = 130.0
    forwardPrice = 125.0 # Locked
    discountRate = 0.05
    dividendRate = 0.02

    ###########################################################################

    expiryDate = valueDate.addMonths(12)
    notional = 100.0

    discountCurve = FinDiscountCurveFlat(valueDate, discountRate)
    dividendCurve = FinDiscountCurveFlat(valueDate, dividendRate)

    equityForward = FinEquityForward(expiryDate,
                                     forwardPrice,
                                     notional,
                                     FinLongShort.LONG)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdPrice = equityForward.forward(valueDate,
                                     stockPrice,
                                     discountCurve, 
                                     dividendCurve)

    fwdValue = equityForward.value(valueDate,
                                   stockPrice,
                                   discountCurve, 
                                   dividendCurve)

#    print(stockPrice, fwdPrice, fwdValue)
    testCases.print(stockPrice, fwdPrice, fwdValue)

###############################################################################


test_FinEquityForward()
testCases.compareTestCases()
