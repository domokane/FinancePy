###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_forward import EquityForward
from financepy.utils.date import Date
from financepy.utils.global_types import FinLongShort
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityForward():

    valuation_date = Date(13, 2, 2018)
    expiry_date = valuation_date.add_months(12)

    stock_price = 130.0
    forward_price = 125.0  # Locked
    discountRate = 0.05
    dividendRate = 0.02

    ###########################################################################

    expiry_date = valuation_date.add_months(12)
    notional = 100.0

    discount_curve = DiscountCurveFlat(valuation_date, discountRate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividendRate)

    equityForward = EquityForward(expiry_date,
                                  forward_price,
                                  notional,
                                  FinLongShort.LONG)

    testCases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwdPrice = equityForward.forward(valuation_date,
                                     stock_price,
                                     discount_curve,
                                     dividend_curve)

    fwdValue = equityForward.value(valuation_date,
                                   stock_price,
                                   discount_curve,
                                   dividend_curve)

#    print(stock_price, fwdPrice, fwdValue)
    testCases.print(stock_price, fwdPrice, fwdValue)

###############################################################################


test_EquityForward()
testCases.compareTestCases()
