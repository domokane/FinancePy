###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date


expiry_date = Date(1, 7, 2015)
call_option = EquityVanillaOption(
    expiry_date, 100.0, OptionTypes.EUROPEAN_CALL
)
call_option_vector = EquityVanillaOption(
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


def test_call_option():
    v = call_option.value(
        value_date, stock_price, discount_curve, dividend_curve, model
    )
    call_option_vector.value(
        value_date, stock_price, discount_curve, dividend_curve, model
    ) == [v] * 3
    assert v.round(4) == 9.3021


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


def test_put_option():
    v = put_option.value(
        value_date, stock_price, discount_curve, dividend_curve, model
    )
    assert v.round(4) == 7.3478
