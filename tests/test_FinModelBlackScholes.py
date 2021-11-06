###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

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
