###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.fx.fx_fixed_lookback_option import FXFixedLookbackOption
from financepy.products.fx.fx_float_lookback_option import FXFloatLookbackOption
from financepy.utils.global_types import OptionTypes


valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.3
num_paths = 10000
stock_priceRange = range(90, 110, 5)
num_steps_per_year = 252

domesticRate = 0.05
domestic_curve = DiscountCurveFlat(valuation_date, domesticRate)

foreignRate = 0.02
foreign_curve = DiscountCurveFlat(valuation_date, foreignRate)


def test_european_call():
    option_type = OptionTypes.EUROPEAN_CALL
    option = FXFloatLookbackOption(expiry_date, option_type)
    stockMin = stock_price - 10
    value = option.value(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin)

    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 23.7455
    assert round(value_mc, 4) == 21.3861

    k = 100.0
    option = FXFixedLookbackOption(expiry_date, option_type, k)
    stockMin = stock_price
    value = option.value(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin)

    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 26.8608
    assert round(value_mc, 4) == 25.6946


def test_european_put():
    option_type = OptionTypes.EUROPEAN_PUT
    option = FXFloatLookbackOption(expiry_date, option_type)
    stockMax = stock_price + 10
    value = option.value(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMax)

    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMax,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 25.2429
    assert round(value_mc, 4) == 25.4885

    k = 105.0
    option = FXFixedLookbackOption(expiry_date, option_type, k)
    stockMin = stock_price - 10.0
    value = option.value(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin)

    value_mc = option.value_mc(
        valuation_date,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stockMin,
        num_paths,
        num_steps_per_year)

    assert round(value, 4) == 25.6047
    assert round(value_mc, 4) == 25.1007
