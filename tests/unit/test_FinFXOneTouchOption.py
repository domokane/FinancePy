# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.global_types import TouchOptionTypes
from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption
from financepy.products.fx.fx_one_touch_option import FXOneTouchOption

value_dt = Date(1, 1, 2016)
expiry_dt = Date(2, 7, 2016)
domestic_curve = FlatDiscountCurve(value_dt, 0.10)
foreign_curve = FlatDiscountCurve(value_dt, 0.03)
model = BlackScholes(0.20)
barrier_rate = 1.0
payment_size = 1.5


########################################################################################


def test_matches_equity_one_touch():

    # An FX rate under Black-Scholes is a stock paying the foreign rate as its dividend yield
    for opt_type in TouchOptionTypes:
        spot_fx_rate = 1.05 if "DOWN" in opt_type.name else 0.95
        fx_option = FXOneTouchOption(expiry_dt, opt_type, barrier_rate, payment_size)
        equity_option = EquityOneTouchOption(expiry_dt, opt_type, barrier_rate, payment_size)
        v_fx = fx_option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)
        v_equity = equity_option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)
        assert np.isclose(v_fx, v_equity, rtol=1e-6), (opt_type, v_fx, v_equity)


########################################################################################


def test_low_volatility():

    up_curves = (FlatDiscountCurve(value_dt, 0.05), FlatDiscountCurve(value_dt, 0.02))
    down_curves = (FlatDiscountCurve(value_dt, 0.02), FlatDiscountCurve(value_dt, 0.05))

    # Expected values are Haug's formulas evaluated in 50-digit arithmetic
    cases = [
        (TouchOptionTypes.UP_AND_IN_CASH_AT_HIT, 1.01, up_curves, 1.4753293742197326),
        (TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY, 1.01, up_curves, 1.4628646684008793),
        (TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING, 1.01, up_curves, 0.0),
        (TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT, 1.01, up_curves, 0.9933884453079533),
        (TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY, 1.01, up_curves, 0.9900227094677558),
        (TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING, 1.01, up_curves, 0.0),
        (TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT, 0.99, down_curves, 1.4899833281515131),
        (TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY, 0.99, down_curves, 1.4850340642016338),
        (TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING, 0.99, down_curves, 0.0),
        (TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT, 0.99, down_curves, 0.9833889965799987),
        (TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY, 0.99, down_curves, 0.9752431122672528),
        (TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING, 0.99, down_curves, 0.0),
    ]

    for opt_type, barrier, curves, expected in cases:
        option = FXOneTouchOption(expiry_dt, opt_type, barrier, payment_size)
        v = option.value(value_dt, 1.0, curves[0], curves[1], BlackScholes(0.0005))
        assert np.isclose(v, expected, rtol=1e-10, atol=1e-10), (opt_type, v, expected)
