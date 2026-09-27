# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.utils.date import Date
import numpy as np

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

    discount_curve = FlatDiscountCurve(value_dt, r)

    use_forward = False

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
    assert round(k1, 4) == 0.0447

    k2 = vol_swap.fair_strike_approx(value_dt, s, strikes, vols)
    assert round(k2, 4) == 0.0424


########################################################################################


def test_fair_strike_with_more_puts_than_positive_strikes():
    """When the requested put strikes would go below zero the replication keeps
    only the positive ones; the weights and the pricing loop must use that
    reduced count instead of raising IndexError. Under a flat volatility the
    fair variance strike approaches the squared volatility."""
    import numpy as np

    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
    from financepy.market.volatility.equity_vol_curve import EquityVolCurve
    from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
    from financepy.utils.date import Date

    value_dt = Date(1, 1, 2026)
    maturity_dt = Date(1, 1, 2027)
    t_exp = (maturity_dt - value_dt) / 365.0
    r, q, sigma, stock_price = 0.05, 0.02, 0.30, 100.0
    discount_curve = FlatDiscountCurve(value_dt, r)
    dividend_curve = FlatDiscountCurve(value_dt, q)
    strikes = np.linspace(20.0, 300.0, 281)
    vol_curve = EquityVolCurve(strikes, np.full(len(strikes), sigma), stock_price, t_exp, r, q)
    swap = EquityVarianceSwap(value_dt, maturity_dt, sigma**2)

    # 140 puts with unit spacing from a forward near 103 reach below zero
    fair_var = swap.fair_strike(
        value_dt, stock_price, dividend_curve, vol_curve, 140, 140, 1.0, discount_curve
    )
    assert swap.num_put_options < 140
    assert len(swap.put_wts) == swap.num_put_options
    assert abs(fair_var - sigma**2) / sigma**2 < 0.02

    # A request that fits is unchanged
    fair_var = swap.fair_strike(
        value_dt, stock_price, dividend_curve, vol_curve, 30, 30, 3.0, discount_curve
    )
    assert swap.num_put_options == 30
    assert abs(fair_var - sigma**2) / sigma**2 < 0.02


def test_fair_strike_with_dividend_yield_and_zero_strike_grid():
    """Under a flat volatility the fair variance strike equals the squared
    volatility whatever the rate and dividend yield. The log-contract drift is
    r - q and the option portfolio compounds at r; a put grid that reaches a zero
    strike must stop before it."""
    import numpy as np

    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
    from financepy.market.volatility.equity_vol_curve import EquityVolCurve
    from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
    from financepy.utils.date import Date

    value_dt = Date(1, 1, 2026)
    maturity_dt = Date(1, 1, 2027)
    t_exp = (maturity_dt - value_dt) / 365.0
    sigma, stock_price = 0.30, 100.0
    strikes = np.linspace(20.0, 300.0, 281)
    swap = EquityVarianceSwap(value_dt, maturity_dt, sigma**2)
    for r, q in [(0.05, 0.02), (0.0, 0.03), (0.05, 0.05), (0.0, 0.0)]:
        discount_curve = FlatDiscountCurve(value_dt, r)
        dividend_curve = FlatDiscountCurve(value_dt, q)
        vol_curve = EquityVolCurve(
            strikes, np.full(len(strikes), sigma), stock_price, t_exp, r, q
        )
        # spacing 0.5 from a forward of 100 reaches a zero strike when r == q
        fair_var = swap.fair_strike(
            value_dt, stock_price, dividend_curve, vol_curve, 200, 200, 0.5, discount_curve
        )
        assert np.isfinite(fair_var), (r, q)
        assert abs(fair_var - sigma**2) / sigma**2 < 0.01, (r, q, fair_var)
