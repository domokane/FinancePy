# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_binomial_tree import (
    EquityTreePayoffTypes,
)
from financepy.products.equity.equity_binomial_tree import (
    EquityTreeExerciseTypes,
)
from financepy.products.equity.equity_binomial_tree import EquityBinomialTree

stock_price = 50.0
risk_free_rate = 0.06
dividend_yield = 0.04
volatility = 0.40

value_dt = Date(1, 1, 2016)
expiry_dt = Date(1, 1, 2017)

model = BlackScholes(volatility)
discount_curve = DiscountCurveFlat(value_dt, risk_free_rate)
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

num_steps = 100

strike_price = 50.0

tree = EquityBinomialTree()

########################################################################################


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
        value_dt,
        payoff,
        expiry_dt,
        payoff,
        exercise,
        params,
    )

    assert [round(x, 3) for x in value] == [7.033, -0.384, 0.019, -2.873]


########################################################################################


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
        value_dt,
        payoff,
        expiry_dt,
        payoff,
        exercise,
        params,
    )

    assert [round(x, 3) for x in value] == [7.219, -0.400, 0.020, -3.128]


########################################################################################


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
        value_dt,
        payoff,
        expiry_dt,
        payoff,
        exercise,
        params,
    )

    assert [round(x, 3) for x in value] == [8.073, 0.577, 0.019, -3.859]


########################################################################################


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
        value_dt,
        payoff,
        expiry_dt,
        payoff,
        exercise,
        params,
    )

    assert [round(x, 3) for x in value] == [8.091, 0.580, 0.019, -3.909]
