###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_american_option import EquityAmericanOption


valuation_date = Date(1, 1, 2016)
expiry_date = Date(1, 1, 2017)
stock_price = 50.0
interest_rate = 0.06
dividend_yield = 0.04
volatility = 0.40
strike_price = 50.0
num_steps = 100

discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

model = BlackScholes(volatility,
                     BlackScholesTypes.CRR_TREE,
                     num_steps)


def test_european_put():
    put_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_PUT)

    value = put_option.value(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)

    assert round(value, 4) == 7.0833


def test_american_put():
    put_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.AMERICAN_PUT)

    value = put_option.value(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)

    assert round(value, 4) == 7.2583


def test_european_call():
    call_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_CALL)

    value = call_option.value(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)

    assert round(value, 4) == 8.0345


def test_american_call():
    call_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.AMERICAN_CALL)

    value = call_option.value(valuation_date, stock_price, discount_curve,
                              dividend_curve, model)

    assert round(value, 4) == 8.0556
