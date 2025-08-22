# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.utils.date import Date
import numpy as np

########################################################################################


def vol_skew(K, atm_vol, atm_k, skew):

    v = atm_vol + skew * (K - atm_k)
    return v

########################################################################################


def test_equity_variance_swap():

    start_dt = Date(20, 3, 2018)
    tenor = "3M"
    strike = 0.3 * 0.3

    vol_swap = EquityVarianceSwap(start_dt, tenor, strike)

    value_dt = Date(20, 3, 2018)
    stock_price = 100.0
    dividend_yield = 0.0
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    maturity_dt = start_dt.add_months(3)

    atm_vol = 0.20
    atm_k = 100.0
    skew = -0.02 / 5.0  # defined as dsigma/dk
    strikes = np.linspace(50.0, 135.0, 18)
    vols = vol_skew(strikes, atm_vol, atm_k, skew)
    vol_curve = EquityVolCurve(value_dt, maturity_dt, strikes, vols)

    strike_spacing = 5.0
    num_call_options = 10
    num_put_options = 10
    r = 0.05

    discount_curve = DiscountCurveFlat(value_dt, r)

    use_forward = False

    k1 = vol_swap.fair_strike(
        value_dt,
        stock_price,
        dividend_curve,
        vol_curve,
        num_call_options,
        num_put_options,
        strike_spacing,
        discount_curve,
        use_forward,
    )
    assert round(k1, 4) == 0.0447

    k2 = vol_swap.fair_strike_approx(value_dt, stock_price, strikes, vols)
    assert round(k2, 4) == 0.0424
