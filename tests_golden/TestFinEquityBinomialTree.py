###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from financepy.products.equity.equity_binomial_tree import EquityBinomialTree
from financepy.products.equity.equity_binomial_tree import EquityTreeExerciseTypes
from financepy.products.equity.equity_binomial_tree import EquityTreePayoffTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.global_types import FinOptionTypes
from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinBinomialTree():

    stock_price = 50.0
    riskFreeRate = 0.06
    dividendYield = 0.04
    volatility = 0.40

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, riskFreeRate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividendYield)

    num_stepsList = [100, 500, 1000, 2000, 5000]

    strike_price = 50.0

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = EquityVanillaOption(
        expiry_date,
        strike_price,
        FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
    delta = putOption.delta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    gamma = putOption.gamma(valuation_date, stock_price, discount_curve, dividend_curve, model)
    theta = putOption.theta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strike_price])

    testCases.header("NumSteps", "Results", "TIME")

    for num_steps in num_stepsList:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
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
        end = time.time()
        duration = end - start
        testCases.print(num_steps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([-1, strike_price])

    testCases.header("NumSteps", "Results", "TIME")

    for num_steps in num_stepsList:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
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
        end = time.time()
        duration = end - start
        testCases.print(num_steps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = EquityVanillaOption(
        expiry_date,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
    delta = callOption.delta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    gamma = callOption.gamma(valuation_date, stock_price, discount_curve, dividend_curve, model)
    theta = callOption.theta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strike_price])

    testCases.header("NumSteps", "Results", "TIME")
    for num_steps in num_stepsList:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
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

        end = time.time()
        duration = end - start
        testCases.print(num_steps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strike_price])

    testCases.header("NumSteps", "Results", "TIME")
    for num_steps in num_stepsList:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
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

        end = time.time()
        duration = end - start
        testCases.print(num_steps, results, duration)

###############################################################################

test_FinBinomialTree()
testCases.compareTestCases()
