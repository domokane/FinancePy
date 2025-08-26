# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_digital_option import FinDigitalOptionTypes
from financepy.products.equity.equity_digital_option import EquityDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date


underlying_type = FinDigitalOptionTypes.CASH_OR_NOTHING

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01
discount_curve = DiscountCurveFlat(value_dt, interest_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

model = BlackScholes(volatility)

num_paths = 40000

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
