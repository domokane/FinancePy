# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.utils.date import Date
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def vol_skew(k, atm_vol, atm_k, skew):

    v = atm_vol + skew * (k - atm_k)
    return v


########################################################################################


def test_equity_variance_swap():

    start_dt = Date(20, 3, 2018)
    tenor = "3M"
    strike = 0.3 * 0.3

    vol_swap = EquityVarianceSwap(start_dt, tenor, strike)

    value_dt = Date(20, 3, 2018)
    s = 100.0
    r = 0.05
    q = 0.0
    dividend_curve = FlatDiscountCurve(value_dt, q)
    discount_curve = FlatDiscountCurve(value_dt, r)

    t = 0.25
    atm_vol = 0.20
    atm_k = 100.0
    skew = -0.02 / 5.0  # defined as dsigma/dk
    strikes = np.linspace(50.0, 135.0, 18)
    vols = vol_skew(strikes, atm_vol, atm_k, skew)
    vol_curve = EquityVolCurve(strikes, vols, s, t, r, q)

    strike_spacing = 5.0
    num_call_options = 10
    num_put_options = 10

    use_forward = False

    test_cases.header("LABEL", "VALUE")

    k1 = vol_swap.fair_strike(
        value_dt,
        s,
        dividend_curve,
        vol_curve,
        num_call_options,
        num_put_options,
        strike_spacing,
        discount_curve,
        use_forward,
    )

    test_cases.print("REPLICATION VARIANCE:", k1)

    k2 = vol_swap.fair_strike_approx(value_dt, s, strikes, vols)
    test_cases.print("DERMAN SKEW APPROX for K:", k2)


########################################################################################

test_equity_variance_swap()
test_cases.compare_test_cases()
