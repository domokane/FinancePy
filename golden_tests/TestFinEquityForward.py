###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

from financepy.products.equity.equity_forward import EquityForward
from financepy.utils.date import Date
from financepy.utils.global_types import FinLongShort
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityForward():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    forward_price = 125.0  # Locked
    discount_rate = 0.05
    dividend_rate = 0.02

    ###########################################################################

    expiry_dt = value_dt.add_months(12)
    notional = 100.0

    discount_curve = DiscountCurveFlat(value_dt, discount_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_rate)

    equityForward = EquityForward(
        expiry_dt, forward_price, notional, FinLongShort.LONG
    )

    test_cases.header("SPOT FX", "FX FWD", "VALUE_BS")

    fwd_price = equityForward.forward(
        value_dt, stock_price, discount_curve, dividend_curve
    )

    fwd_value = equityForward.value(
        value_dt, stock_price, discount_curve, dividend_curve
    )

    #    print(stock_price, fwd_price, fwd_value)
    test_cases.print(stock_price, fwd_price, fwd_value)


###############################################################################


test_EquityForward()
test_cases.compareTestCases()
