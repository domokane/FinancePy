###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.fx.fx_fixed_lookback_option import (
    FXFixedLookbackOption,
)
from financepy.products.fx.fx_float_lookback_option import (
    FXFloatLookbackOption,
)
from financepy.utils.global_types import OptionTypes


value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.3
num_paths = 10000
stock_priceRange = range(90, 110, 5)
num_steps_per_year = 252

domesticRate = 0.05
domestic_curve = DiscountCurveFlat(value_dt, domesticRate)

foreignRate = 0.02
foreign_curve = DiscountCurveFlat(value_dt, foreignRate)


def test_european_call():
    option_type = OptionTypes.EUROPEAN_CALL
    option = FXFloatLookbackOption(expiry_dt, option_type)
    stock_min = stock_price - 10
    value = option.value(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
        num_paths,
        num_steps_per_year,
    )

    assert round(value, 4) == 23.7455
    assert round(value_mc, 4) == 23.2553

    k = 100.0
    option = FXFixedLookbackOption(expiry_dt, option_type, k)
    stock_min = stock_price
    value = option.value(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
        num_paths,
        num_steps_per_year,
    )

    assert round(value, 4) == 26.8608
    assert round(value_mc, 4) == 25.7191


def test_european_put():
    option_type = OptionTypes.EUROPEAN_PUT
    option = FXFloatLookbackOption(expiry_dt, option_type)
    stock_max = stock_price + 10
    value = option.value(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_max,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_max,
        num_paths,
        num_steps_per_year,
    )

    assert round(value, 4) == 25.2429
    assert round(value_mc, 4) == 24.3641

    k = 105.0
    option = FXFixedLookbackOption(expiry_dt, option_type, k)
    stock_min = stock_price - 10.0
    value = option.value(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        domestic_curve,
        foreign_curve,
        volatility,
        stock_min,
        num_paths,
        num_steps_per_year,
    )

    assert round(value, 4) == 25.6047
    assert round(value_mc, 4) == 25.0960
