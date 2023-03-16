###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from pytest import approx

from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.models.black_scholes import BlackScholesTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import OptionTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date


valuation_date = Date(8, 5, 2015)
expiry_date = Date(15, 1, 2016)

strike_price = 130.0
stock_price = 127.62
volatility = 0.20
interest_rate = 0.001
dividend_yield = 0.0163

option_type = OptionTypes.AMERICAN_CALL
euOptionType = OptionTypes.EUROPEAN_CALL

amOption = EquityAmericanOption(expiry_date, strike_price,
                                option_type)

ameuOption = EquityAmericanOption(expiry_date, strike_price,
                                  euOptionType)

euOption = EquityVanillaOption(expiry_date, strike_price,
                               euOptionType)

discount_curve = DiscountCurveFlat(valuation_date, interest_rate,
                                   FrequencyTypes.CONTINUOUS,
                                   DayCountTypes.ACT_365F)

dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield,
                                   FrequencyTypes.CONTINUOUS,
                                   DayCountTypes.ACT_365F)

num_steps_per_year = 400

modelTree = BlackScholes(volatility,
                         BlackScholesTypes.CRR_TREE,
                         num_steps_per_year)


def test_black_scholes():
    v = amOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, modelTree)
    assert round(v, 4) == 6.8391

    modelApprox = BlackScholes(volatility,
                               BlackScholesTypes.BARONE_ADESI)

    v = amOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, modelApprox)

    assert round(v, 4) == 6.8277

    v = ameuOption.value(valuation_date, stock_price, discount_curve,
                         dividend_curve, modelTree)

    assert round(v, 4) == 6.7510

    v = euOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, modelTree)

    assert round(v, 4) == 6.7493


def test_bjerksund_stensland():
    # Valuation of American call option as in Bjerksund and Sensland's paper published in 1993.
    # See Table 2 in https://www.sciencedirect.com/science/article/abs/pii/095652219390009H

    # valuation_date and exipry_date are set so that time to maturity becomes 0.25
    valuation_date = Date(8, 5, 2015)
    expiry_date = Date(7, 8, 2015, hh=6)
    interest_rate = 0.08
    volatility = 0.40
    borrow_rate = 0.04
    strike_price = 100
    stock_prices = [80.0, 90.0, 100.0, 110.0, 120.0]

    # model setting
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate,
                                       FrequencyTypes.CONTINUOUS,
                                       DayCountTypes.ACT_365F)
    borrow_curve = DiscountCurveFlat(valuation_date, borrow_rate,
                                     FrequencyTypes.CONTINUOUS,
                                     DayCountTypes.ACT_365F)
    model = BlackScholes(volatility, BlackScholesTypes.Bjerksund_Stensland)

    # american call case
    amCallOption = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.AMERICAN_CALL)
    values = []
    for stock_price in stock_prices:
        value = amCallOption.value(valuation_date, stock_price, discount_curve,
                                   borrow_curve, model)
        values.append(round(value, 2))
    assert values == [1.29, 3.82, 8.35, 14.80, 22.71]

    # american put case
    amPutOption = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.AMERICAN_PUT)
    values = []
    for stock_price in stock_prices:
        value = amPutOption.value(valuation_date, stock_price, discount_curve,
                                  borrow_curve, model)
        values.append(round(value, 2))
    assert values == [20.53, 12.91, 7.42, 3.93, 1.93]


def test_black_scholes_finite_difference():
    """
    Assert finite difference model matches tree model to at least 1 dp
    """
    params = {
        'num_samples': 200,
        'theta': 0.5
    }
    model = BlackScholes(volatility,
                         implementationType=BlackScholesTypes.FINITE_DIFFERENCE,
                         params=params)

    v = amOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

    assert v == approx(6.8391, 1e-1)

    v = ameuOption.value(valuation_date, stock_price, discount_curve,
                         dividend_curve, model)

    assert v == approx(6.7510, 1e-1)

    v = euOption.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

    assert v == approx(6.7493, 1e-1)
