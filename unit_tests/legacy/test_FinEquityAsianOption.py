# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.equity.equity_asian_option import (
    AsianOptionValuationMethods,
)
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.utils.global_types import OptionTypes


value_dt = Date(1, 1, 2014)
start_averaging_dt = Date(1, 6, 2014)
expiry_dt = Date(1, 1, 2015)
stock_price = 100.0
volatility = 0.20
interest_rate = 0.30
dividend_yield = 0.10
num_observations = 120  # daily as we have a half year
accrued_avg = None
k = 100
seed = 1976
num_paths = 5000

model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(value_dt, interest_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

asian_option = EquityAsianOption(
    start_averaging_dt,
    expiry_dt,
    k,
    OptionTypes.EUROPEAN_CALL,
    num_observations,
)

########################################################################################


def test_geometric():

    value_geometric = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationMethods.GEOMETRIC,
        accrued_avg,
    )

    assert round(value_geometric, 4) == 12.3380


########################################################################################


def test_turnbull_wakeman():

    value_turnbull_wakeman = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationMethods.TURNBULL_WAKEMAN,
        accrued_avg,
    )

    assert round(value_turnbull_wakeman, 4) == 12.5381


########################################################################################


def test_curran():

    value_curran = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationMethods.CURRAN,
        accrued_avg,
    )

    assert round(value_curran, 4) == 12.5368


########################################################################################


def test_mc():

    value_mc = asian_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
        seed,
        accrued_avg,
    )

    assert round(value_mc, 4) == 12.5722
