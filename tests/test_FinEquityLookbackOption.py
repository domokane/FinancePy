###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_fixed_lookback_option import EquityFixedLookbackOption
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption


valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.3
interest_rate = 0.05
dividend_yield = 0.01
num_paths = 10000
stock_priceRange = range(90, 110, 10)
num_steps_per_year = 252

discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)


def test_european_call():
    option_type = OptionTypes.EUROPEAN_CALL
    k = 100.0
    option = EquityFixedLookbackOption(expiry_date, option_type, k)

    stockMax = stock_price + 10.0
    value = option.value(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        stockMax)
    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        stockMax,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 28.7477
    assert round(value_mc, 4) == 27.8592


def test_european_put():
    option_type = OptionTypes.EUROPEAN_PUT
    k = 100.0
    option = EquityFixedLookbackOption(expiry_date, option_type, k)

    stockMin = stock_price - 10
    value = option.value(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        stockMin)
    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        stockMin,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 20.5366
    assert round(value_mc, 4) == 20.0366
