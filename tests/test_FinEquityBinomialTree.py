###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_binomial_tree import EquityTreePayoffTypes
from financepy.products.equity.equity_binomial_tree import EquityTreeExerciseTypes
from financepy.products.equity.equity_binomial_tree import EquityBinomialTree
import numpy as np

stock_price = 50.0
risk_free_rate = 0.06
dividend_yield = 0.04
volatility = 0.40

valuation_date = Date(1, 1, 2016)
expiry_date = Date(1, 1, 2017)

model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(valuation_date, risk_free_rate)
dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

num_steps = 100

strike_price = 50.0

tree = EquityBinomialTree()


def test_european_put():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strike_price])

    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)

    assert [round(x, 4) for x in value] == [7.1050, -0.3865, 0.0187, -2.9453]


def test_american_put():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([-1, strike_price])

    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)

    assert [round(x, 4) for x in value] == [7.2753, -0.4008, 0.0200, -3.1803]


def test_european_call():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strike_price])

    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)

    assert [round(x, 4) for x in value] == [8.0175, 0.5747, 0.0187, -3.8111]


def test_american_call():
    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strike_price])

    value = tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        num_steps,
        valuation_date,
        payoff,
        expiry_date,
        payoff,
        exercise,
        params)

    assert [round(x, 4) for x in value] == [8.0399, 0.5775, 0.0189, -3.8676]
