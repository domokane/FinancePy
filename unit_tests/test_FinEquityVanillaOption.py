# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import math

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date


expiry_date = Date(1, 7, 2015)
call_option = EquityVanillaOption(expiry_date, 100.0, OptionTypes.EUROPEAN_CALL)
call_optionormcdf_vector = EquityVanillaOption(
    [expiry_date] * 3, 100.0, OptionTypes.EUROPEAN_CALL
)
put_option = EquityVanillaOption(expiry_date, 100.0, OptionTypes.EUROPEAN_PUT)

value_date = Date(1, 1, 2015)
stock_price = 100
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01
model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(value_date, interest_rate)
dividend_curve = DiscountCurveFlat(value_date, dividend_yield)

########################################################################################


def test_call_option():

    v = call_option.value(
        value_date, stock_price, discount_curve, dividend_curve, model
    )
    call_optionormcdf_vector.value(
        value_date, stock_price, discount_curve, dividend_curve, model
    ) == [v] * 3
    assert v.round(4) == 9.3021


########################################################################################


def test_greeks():

    delta = call_option.delta(
        value_date, stock_price, discount_curve, dividend_curve, model
    )

    vega = call_option.vega(
        value_date, stock_price, discount_curve, dividend_curve, model
    )

    theta = call_option.theta(
        value_date, stock_price, discount_curve, dividend_curve, model
    )

    rho = call_option.rho(
        value_date, stock_price, discount_curve, dividend_curve, model
    )

    assert [round(x, 4) for x in (delta, vega, theta, rho)] == [
        0.5762,
        27.4034,
        -10.1289,
        23.9608,
    ]


########################################################################################


def test_put_option():

    v = put_option.value(value_date, stock_price, discount_curve, dividend_curve, model)
    assert v.round(4) == 7.3478


########################################################################################


def closed_form_vanna(s, k, t, r, q, v):
    # Black-Scholes vanna = d(vega)/d(spot) = -exp(-q t) n(d1) d2 / v
    d1 = (math.log(s / k) + (r - q + 0.5 * v * v) * t) / (v * math.sqrt(t))
    d2 = d1 - v * math.sqrt(t)
    n_d1 = math.exp(-0.5 * d1 * d1) / math.sqrt(2.0 * math.pi)
    return -math.exp(-q * t) * n_d1 * d2 / v


def fd_vanna(option, s):
    # d(delta)/d(vol) by central difference of the library's own delta
    h = 1e-4
    up = option.delta(
        value_date, s, discount_curve, dividend_curve, BlackScholes(volatility + h)
    )
    dn = option.delta(
        value_date, s, discount_curve, dividend_curve, BlackScholes(volatility - h)
    )
    return (up - dn) / (2 * h)


def test_vanna():

    t = (expiry_date - value_date) / 365.0

    # OTM / ATM / ITM: vanna is positive below the strike and negative above it
    for s, expected in ((80, 0.9791), (100, 0.0152), (110, -0.4705)):

        call = EquityVanillaOption(expiry_date, 100.0, OptionTypes.EUROPEAN_CALL)
        put = EquityVanillaOption(expiry_date, 100.0, OptionTypes.EUROPEAN_PUT)

        v_call = float(call.vanna(value_date, s, discount_curve, dividend_curve, model))
        v_put = float(put.vanna(value_date, s, discount_curve, dividend_curve, model))

        ref = closed_form_vanna(s, 100.0, t, interest_rate, dividend_yield, volatility)

        assert round(v_call, 4) == expected
        assert abs(v_call - ref) < 1e-6
        assert abs(v_call - fd_vanna(call, s)) < 1e-3
        assert abs(v_call - v_put) < 1e-9  # vanna does not depend on option type
