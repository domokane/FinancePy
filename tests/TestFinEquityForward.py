###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_forward import EquityForward
from financepy.utils.date import Date
from financepy.utils.global_types import FinLongShort
from financepy.market.discount.curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityForward():

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
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendRate)

    equityForward = EquityForward(expiryDate,
                                     forwardPrice,
                                     notional,
                                     FinLongShort.LONG)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdPrice = equityForward.forward(valueDate,
                                     stockPrice,
                                     discountCurve, 
                                     dividend_curve)

    fwdValue = equityForward.value(valueDate,
                                   stockPrice,
                                   discountCurve, 
                                   dividend_curve)

#    print(stockPrice, fwdPrice, fwdValue)
    testCases.print(stockPrice, fwdPrice, fwdValue)

###############################################################################


test_EquityForward()
testCases.compareTestCases()
