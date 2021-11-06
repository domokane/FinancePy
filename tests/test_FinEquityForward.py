###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import FinLongShort
from financepy.utils.date import Date
from financepy.products.equity.equity_forward import EquityForward


def test_equity_forward():

    valuation_date = Date(13, 2, 2018)
    expiry_date = valuation_date.add_months(12)

    stock_price = 130.0
    forward_price = 125.0  # Locked
    discountRate = 0.05
    dividendRate = 0.02

    expiry_date = valuation_date.add_months(12)
    notional = 100.0

    discount_curve = DiscountCurveFlat(valuation_date, discountRate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividendRate)

    equityForward = EquityForward(expiry_date,
                                  forward_price,
                                  notional,
                                  FinLongShort.LONG)

    fwdPrice = equityForward.forward(valuation_date,
                                     stock_price,
                                     discount_curve,
                                     dividend_curve)

    fwdValue = equityForward.value(valuation_date,
                                   stock_price,
                                   discount_curve,
                                   dividend_curve)

    assert round(fwdPrice, 4) == 133.9591
    assert round(fwdValue, 4) == 852.2149
