# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.products.equity.equity_compound_option import EquityCompoundOption
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

    assert round(value, 4) == 4.6039
    assert round(values[0], 4) == 4.5982

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

    assert round(value, 4) == 6.7176
    assert round(values[0], 4) == 6.7053

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

    assert round(value, 4) == 2.0165
    assert round(values[0], 4) == 2.0228

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

    assert round(value, 4) == 1.0873
    assert round(values[0], 4) == 1.0869


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

    assert round(value, 4) == 4.5982
    assert round(values[0], 4) == 4.5982

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

    assert round(value, 4) == 12.8406
    assert round(values[0], 4) == 12.8406

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

    assert round(value, 4) == 4.6652
    assert round(values[0], 4) == 4.6652

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

    assert round(value, 4) == 4.2982
    assert round(values[0], 4) == 4.2982


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

    assert round(delta, 4) == -0.7191
    assert round(vega, 4) == 0.3758
    assert round(theta, 4) == 0.7478
