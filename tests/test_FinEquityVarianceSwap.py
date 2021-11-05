###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.utils.date import Date
import numpy as np


def vol_skew(K, atm_vol, atmK, skew):
    v = atm_vol + skew * (K-atmK)
    return v


def test_equity_variance_swap():
    start_date = Date(20, 3, 2018)
    tenor = "3M"
    strike = 0.3*0.3

    volSwap = EquityVarianceSwap(start_date, tenor, strike)

    valuation_date = Date(20, 3, 2018)
    stock_price = 100.0
    dividend_yield = 0.0
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    maturity_date = start_date.add_months(3)

    atm_vol = 0.20
    atmK = 100.0
    skew = -0.02/5.0  # defined as dsigma/dK
    strikes = np.linspace(50.0, 135.0, 18)
    vols = vol_skew(strikes, atm_vol, atmK, skew)
    volCurve = EquityVolCurve(valuation_date, maturity_date, strikes, vols)

    strike_spacing = 5.0
    num_call_options = 10
    num_put_options = 10
    r = 0.05

    discount_curve = DiscountCurveFlat(valuation_date, r)

    use_forward = False

    k1 = volSwap.fair_strike(valuation_date, stock_price, dividend_curve,
                             volCurve, num_call_options, num_put_options,
                             strike_spacing, discount_curve, use_forward)
    assert round(k1, 4) == 0.0447

    k2 = volSwap.fair_strike_approx(valuation_date, stock_price, strikes, vols)
    assert round(k2, 4) == 0.0424
