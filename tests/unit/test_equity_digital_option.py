# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.global_types import OptionTypes
from financepy.utils.global_types import DigitalOptionTypes
from financepy.products.equity.equity_digital_option import EquityDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date


underlying_type = DigitalOptionTypes.CASH_OR_NOTHING

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

model = BlackScholes(volatility)

num_paths = 40000


@pytest.mark.parametrize("call_put", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT])
@pytest.mark.parametrize("digital_type", [DigitalOptionTypes.CASH_OR_NOTHING, DigitalOptionTypes.ASSET_OR_NOTHING])
@pytest.mark.parametrize("spot", [90.0, 100.0, 110.0])
def test_expiry_analytic_and_mc_match_strict_event_payoff(call_put, digital_type, spot):
    """At expiry neither path simulates a future crossing of the barrier."""
    curve = FlatDiscountCurve(expiry_dt, 0.05)
    option = EquityDigitalOption(expiry_dt, 100.0, call_put, digital_type)
    event = spot > 100.0 if call_put == OptionTypes.EUROPEAN_CALL else spot < 100.0
    expected = float(event)
    if digital_type == DigitalOptionTypes.ASSET_OR_NOTHING:
        expected *= spot
    assert option.value(expiry_dt, spot, curve, curve, model) == expected
    assert option.value_mc(expiry_dt, spot, curve, curve, model) == expected


@pytest.mark.parametrize("call_put", [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT])
@pytest.mark.parametrize("digital_type", [DigitalOptionTypes.CASH_OR_NOTHING, DigitalOptionTypes.ASSET_OR_NOTHING])
def test_expiry_analytic_preserves_vector_payoffs(call_put, digital_type):
    """Analytic vector inputs retain one deterministic payoff per stock."""
    curve = FlatDiscountCurve(expiry_dt, 0.05)
    option = EquityDigitalOption(expiry_dt, 100.0, call_put, digital_type)
    spots = np.array([90.0, 100.0, 110.0])
    expected = np.array([0.0, 0.0, 1.0]) if call_put == OptionTypes.EUROPEAN_CALL else np.array([1.0, 0.0, 0.0])
    if digital_type == DigitalOptionTypes.ASSET_OR_NOTHING:
        expected *= spots
    np.testing.assert_array_equal(option.value(expiry_dt, spots, curve, curve, model), expected)

########################################################################################


def test_value():

    call_option = EquityDigitalOption(
        expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type
    )
    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    value_mc = call_option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_paths
    )

    assert round(value, 4) == 0.4693
    assert round(value_mc, 4) == 0.4694


########################################################################################


def test_greeks():

    call_option = EquityDigitalOption(
        expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type
    )

    delta = call_option.delta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    vega = call_option.vega(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    theta = call_option.theta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(delta, 4) == 0.0126
    assert round(vega, 4) == -0.0035
    assert round(theta, 4) == 0.0266
