###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

sys.path.append("..")

import numpy as np

from financepy.utils.date import Date
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def volSkew(K, atm_vol, atmK, skew):
    v = atm_vol + skew * (K - atmK)
    return v


###############################################################################


def test_EquityVarianceSwap():

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
    atmK = 100.0
    skew = -0.02 / 5.0  # defined as dsigma/dK
    strikes = np.linspace(50.0, 135.0, 18)
    vols = volSkew(strikes, atm_vol, atmK, skew)
    vol_curve = EquityVolCurve(value_dt, maturity_dt, strikes, vols)

    strike_spacing = 5.0
    num_call_options = 10
    num_put_options = 10
    r = 0.05

    discount_curve = DiscountCurveFlat(value_dt, r)

    use_forward = False

    test_cases.header("LABEL", "VALUE")

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

    test_cases.print("REPLICATION VARIANCE:", k1)

    k2 = vol_swap.fair_strike_approx(value_dt, stock_price, strikes, vols)
    test_cases.print("DERMAN SKEW APPROX for K:", k2)


##########################################################################


test_EquityVarianceSwap()
test_cases.compareTestCases()
