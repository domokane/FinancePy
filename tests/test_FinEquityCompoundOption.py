###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
import sys
sys.path.append("./..")


valuation_date = Date(1, 1, 2015)
expiry_date1 = Date(1, 1, 2017)
expiry_date2 = Date(1, 1, 2018)
k1 = 5.0
k2 = 95.0
stock_price = 85.0
volatility = 0.15
interest_rate = 0.035
dividend_yield = 0.01

model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

num_steps = 200


def test_european():
    stock_price = 85.0

    option_type1 = OptionTypes.EUROPEAN_CALL
    option_type2 = OptionTypes.EUROPEAN_CALL

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 4.6039
    assert round(values[0], 4) == 4.5587

    option_type1 = OptionTypes.EUROPEAN_CALL
    option_type2 = OptionTypes.EUROPEAN_PUT

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 6.7176
    assert round(values[0], 4) == 6.7676

    option_type1 = OptionTypes.EUROPEAN_PUT
    option_type2 = OptionTypes.EUROPEAN_CALL

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 2.0165
    assert round(values[0], 4) == 2.0361

    option_type1 = OptionTypes.EUROPEAN_PUT
    option_type2 = OptionTypes.EUROPEAN_PUT

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 1.0873
    assert round(values[0], 4) == 1.0789


def test_american():
    stock_price = 85.0

    option_type1 = OptionTypes.AMERICAN_CALL
    option_type2 = OptionTypes.AMERICAN_CALL

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 4.5587
    assert round(values[0], 4) == 4.5587

    option_type1 = OptionTypes.AMERICAN_CALL
    option_type2 = OptionTypes.AMERICAN_PUT

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 12.8630
    assert round(values[0], 4) == 12.8630

    option_type1 = OptionTypes.AMERICAN_PUT
    option_type2 = OptionTypes.AMERICAN_CALL

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 4.6697
    assert round(values[0], 4) == 4.6697

    option_type1 = OptionTypes.AMERICAN_PUT
    option_type2 = OptionTypes.AMERICAN_PUT

    cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                      expiry_date2, option_type2, k2)
    value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
    values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                    dividend_curve, model, num_steps)

    assert round(value, 4) == 4.3034
    assert round(values[0], 4) == 4.3034


def test_greeks():
    stock_price = 70
    option_type1 = OptionTypes.EUROPEAN_CALL
    option_type2 = OptionTypes.EUROPEAN_PUT
    cmpdOption = EquityCompoundOption(
        expiry_date1, option_type1, k1,
        expiry_date2, option_type2, k2)

    delta = cmpdOption.delta(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        model)
    vega = cmpdOption.vega(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        model)
    theta = cmpdOption.theta(
        valuation_date,
        stock_price,
        discount_curve,
        dividend_curve,
        model)

    assert round(delta, 4) == -0.7191
    assert round(vega, 4) == 0.3758
    assert round(theta, 4) == 0.7478
