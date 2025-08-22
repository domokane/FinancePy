# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
import numpy as np
from financepy.products.equity.equity_binomial_tree import EquityBinomialTree
from financepy.products.equity.equity_binomial_tree import (
    EquityTreeExerciseTypes,
)
from financepy.products.equity.equity_binomial_tree import (
    EquityTreePayoffTypes,
)
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)



########################################################################################


def test_FinBinomialTree():

    stock_price = 50.0
    risk_free_rate = 0.06
    dividend_yield = 0.04
    volatility = 0.40

    value_dt = Date(1, 1, 2016)
    expiry_dt = Date(1, 1, 2017)

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, risk_free_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    num_steps_list = [100, 500, 1000]

    strike_price = 50.0

    test_cases.banner("================== EUROPEAN PUT =======================")

    put_option = EquityVanillaOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT)
    value = put_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    delta = put_option.delta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    gamma = put_option.gamma(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    theta = put_option.theta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    test_cases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    test_cases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strike_price])

    test_cases.header("num_steps", "Results", "TIME")

    for num_steps in num_steps_list:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
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
        end = time.time()
        duration = end - start
        test_cases.print(num_steps, results, duration)

    test_cases.banner("================== AMERICAN PUT =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([-1, strike_price])

    test_cases.header("num_steps", "Results", "TIME")

    for num_steps in num_steps_list:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
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
        end = time.time()
        duration = end - start
        test_cases.print(num_steps, results, duration)

    test_cases.banner("================== EUROPEAN CALL =======================")

    call_option = EquityVanillaOption(
        expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
    )
    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    delta = call_option.delta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    gamma = call_option.gamma(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    theta = call_option.theta(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )
    test_cases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    test_cases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strike_price])

    test_cases.header("num_steps", "Results", "TIME")
    for num_steps in num_steps_list:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
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

        end = time.time()
        duration = end - start
        test_cases.print(num_steps, results, duration)

    test_cases.banner("================== AMERICAN CALL =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strike_price])

    test_cases.header("num_steps", "Results", "TIME")
    for num_steps in num_steps_list:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
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

        end = time.time()
        duration = end - start
        test_cases.print(num_steps, results, duration)


########################################################################################


test_FinBinomialTree()
test_cases.compare_test_cases()

########################################################################################

