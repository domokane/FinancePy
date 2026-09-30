########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.products.fx.fx_variance_swap import FinFXVarianceSwap
from financepy.utils.date import Date

value_dt = Date(1, 1, 2026)
maturity_dt = Date(1, 1, 2027)
t_exp = (maturity_dt - value_dt) / 365.0
sigma = 0.12
spot_fx_rate = 1.10
strikes = np.linspace(0.3, 2.5, 221)

########################################################################################


def _fair_strike(r_dom, r_for, num_call_options=100, num_put_options=100, spacing=0.005):
    domestic_curve = FlatDiscountCurve(value_dt, r_dom)
    foreign_curve = FlatDiscountCurve(value_dt, r_for)
    vol_curve = EquityVolCurve(
        strikes, np.full(len(strikes), sigma), spot_fx_rate, t_exp, r_dom, r_for
    )
    swap = FinFXVarianceSwap(value_dt, maturity_dt, sigma**2)
    fair_var = swap.fair_strike(
        value_dt,
        spot_fx_rate,
        foreign_curve,
        vol_curve,
        num_call_options,
        num_put_options,
        spacing,
        domestic_curve,
    )
    return swap, fair_var


def test_fair_strike_equals_flat_variance_for_any_rates():
    """Under a flat volatility the fair variance strike is the squared volatility
    whatever the domestic and foreign rates. The foreign rate plays the role of a
    dividend yield: the log-contract drift is r_d - r_f and the option portfolio
    compounds at r_d."""
    for r_dom, r_for in [(0.03, 0.0), (0.03, 0.05), (0.01, 0.04), (0.04, 0.04), (0.0, 0.0)]:
        _, fair_var = _fair_strike(r_dom, r_for)
        assert np.isfinite(fair_var), (r_dom, r_for)
        assert abs(fair_var - sigma**2) / sigma**2 < 0.01, (r_dom, r_for, fair_var)


def test_fair_strike_with_more_puts_than_positive_strikes():
    """A put grid that would reach below zero is truncated to positive strikes and
    the weights and pricing loop use the reduced count."""
    swap, fair_var = _fair_strike(0.03, 0.05, num_put_options=300)
    assert swap.num_put_options < 300
    assert len(swap.put_wts) == swap.num_put_options
    assert np.all(np.asarray(swap.put_strikes) > 0.0)
    assert abs(fair_var - sigma**2) / sigma**2 < 0.01


def test_realised_variance_of_constant_log_returns():
    """A series growing by the same log return every day has a realised variance of
    252 times that squared return."""
    swap = FinFXVarianceSwap(value_dt, maturity_dt, sigma**2)
    daily = 0.01
    closes = spot_fx_rate * np.exp(daily * np.arange(11))
    assert np.isclose(swap.realised_variance(closes), 252.0 * daily**2)
