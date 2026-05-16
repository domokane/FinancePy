# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.products.equity.equity_compound_option import (
    EquityCompoundOption,
)
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date


value_dt = Date(1, 1, 2015)
expiry_dt1 = Date(1, 1, 2017)
expiry_dt2 = Date(1, 1, 2018)
k1 = 5.0
k2 = 95.0
stock_price = 85.0
volatility = 0.15
interest_rate = 0.035
dividend_yield = 0.01

model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(value_dt, interest_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

num_steps = 200

########################################################################################


def test_european():

    stock_price = 85.0

    opt_type1 = OptionTypes.EUROPEAN_CALL
    opt_type2 = OptionTypes.EUROPEAN_CALL

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 4.606
    assert round(values[0], 3) == 4.600

    opt_type1 = OptionTypes.EUROPEAN_CALL
    opt_type2 = OptionTypes.EUROPEAN_PUT

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 6.714
    assert round(values[0], 3) == 6.702

    opt_type1 = OptionTypes.EUROPEAN_PUT
    opt_type2 = OptionTypes.EUROPEAN_CALL

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 2.016
    assert round(values[0], 3) == 2.022

    opt_type1 = OptionTypes.EUROPEAN_PUT
    opt_type2 = OptionTypes.EUROPEAN_PUT

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 1.087
    assert round(values[0], 3) == 1.087


########################################################################################


def test_american():

    stock_price = 85.0

    opt_type1 = OptionTypes.AMERICAN_CALL
    opt_type2 = OptionTypes.AMERICAN_CALL

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 4.600
    assert round(values[0], 3) == 4.600

    opt_type1 = OptionTypes.AMERICAN_CALL
    opt_type2 = OptionTypes.AMERICAN_PUT

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 2) == 12.84
    assert round(values[0], 2) == 12.84

    opt_type1 = OptionTypes.AMERICAN_PUT
    opt_type2 = OptionTypes.AMERICAN_CALL

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 4.665
    assert round(values[0], 3) == 4.665

    opt_type1 = OptionTypes.AMERICAN_PUT
    opt_type2 = OptionTypes.AMERICAN_PUT

    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )
    value = cmpd_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    values = cmpd_option.value_tree(
        value_dt, stock_price, discount_curve, dividend_curve, model, num_steps
    )

    assert round(value, 3) == 4.298
    assert round(values[0], 3) == 4.298


########################################################################################


def test_greeks():

    stock_price = 70
    opt_type1 = OptionTypes.EUROPEAN_CALL
    opt_type2 = OptionTypes.EUROPEAN_PUT
    cmpd_option = EquityCompoundOption(
        expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
    )

    delta = cmpd_option.delta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    vega = cmpd_option.vega(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    theta = cmpd_option.theta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    assert round(delta, 3) == -0.719
    assert round(vega, 3) == 0.376
    assert round(theta, 3) == 0.747
