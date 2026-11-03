# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_types import DoubleBarrierTypes
from financepy.products.fx.fx_double_one_touch_option import FXDoubleOneTouchOption
from financepy.models.bs_double_touch_option import fast_double_no_touch_pricer

value_dt = Date(1, 1, 2016)
expiry_dt = Date(1, 7, 2016)
domestic_curve = FlatDiscountCurve(value_dt, 0.05)
foreign_curve = FlatDiscountCurve(value_dt, 0.02)
model = BlackScholes(0.30)
lower = 85.0
upper = 118.0
payment_size = 1.0

knock_out = FXDoubleOneTouchOption(expiry_dt, DoubleBarrierTypes.KNOCK_OUT, lower, upper, payment_size)
knock_in = FXDoubleOneTouchOption(expiry_dt, DoubleBarrierTypes.KNOCK_IN, lower, upper, payment_size)


def _values(spot):
    v_out = knock_out.value(value_dt, spot, domestic_curve, foreign_curve, model)
    v_in = knock_in.value(value_dt, spot, domestic_curve, foreign_curve, model)
    return v_out, v_in


########################################################################################


def test_no_touch_plus_one_touch_is_discounted_payment():

    df_d = domestic_curve.df(expiry_dt)

    for spot in [86.0, 90.0, 100.0, 110.0, 117.0]:
        v_out, v_in = _values(spot)
        assert 0.0 < v_out < payment_size * df_d
        assert 0.0 < v_in < payment_size * df_d
        assert v_out + v_in == pytest.approx(payment_size * df_d, abs=1e-12)


def test_spot_beyond_a_barrier_has_already_touched():

    # The double no-touch is worthless and the double one-touch pays K at expiry.
    df_d = domestic_curve.df(expiry_dt)

    for spot in [80.0, 125.0]:
        v_out, v_in = _values(spot)
        assert v_out == 0.0
        assert v_in == pytest.approx(payment_size * df_d, abs=1e-12)


def test_spot_on_a_barrier_counts_as_a_touch():

    df_d = domestic_curve.df(expiry_dt)

    for spot in [lower, upper]:
        v_out, v_in = _values(spot)
        assert v_out == 0.0
        assert v_in == pytest.approx(payment_size * df_d, abs=1e-12)


def test_values_are_continuous_at_the_barriers():

    df_d = domestic_curve.df(expiry_dt)

    for spot in [lower * (1.0 + 1e-9), upper * (1.0 - 1e-9)]:
        v_out, v_in = _values(spot)
        assert v_out == pytest.approx(0.0, abs=1e-6)
        assert v_in == pytest.approx(payment_size * df_d, abs=1e-6)


def test_fast_pricer_returns_the_no_touch_value_for_both_types():

    t_exp = 0.5
    r_d = 0.05
    r_f = 0.02
    df_d = np.exp(-r_d * t_exp)

    for opt_type in [DoubleBarrierTypes.KNOCK_IN.value, DoubleBarrierTypes.KNOCK_OUT.value]:
        # already touched
        assert fast_double_no_touch_pricer(80.0, lower, upper, 1.0, t_exp, opt_type, r_d, r_f, 0.3) == 0.0
        assert fast_double_no_touch_pricer(upper, lower, upper, 1.0, t_exp, opt_type, r_d, r_f, 0.3) == 0.0
        # deterministic path that stays inside
        assert fast_double_no_touch_pricer(100.0, lower, upper, 1.0, t_exp, opt_type, r_d, r_f, 0.0) == pytest.approx(df_d)
        # deterministic path that drifts through the upper barrier
        assert fast_double_no_touch_pricer(117.0, lower, upper, 1.0, 20.0, opt_type, r_d, r_f, 0.0) == 0.0
