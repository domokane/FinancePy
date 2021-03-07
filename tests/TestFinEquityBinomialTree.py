###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityBinomialTree import FinEquityBinomialTree
from financepy.products.equity.FinEquityBinomialTree import FinEquityTreeExerciseTypes
from financepy.products.equity.FinEquityBinomialTree import FinEquityTreePayoffTypes
from financepy.products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from financepy.utils.global_types import FinOptionTypes
from financepy.utils.date import Date
from financepy.models.black_scholes import FinModelBlackScholes
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

    model = FinModelBlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, riskFreeRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    num_stepsList = [100, 500, 1000, 2000, 5000]

    strikePrice = 50.0

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = FinEquityVanillaOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
    delta = putOption.delta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    gamma = putOption.gamma(valuation_date, stock_price, discount_curve, dividendCurve, model)
    theta = putOption.theta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = FinEquityTreePayoffTypes.VANILLA_OPTION
    exercise = FinEquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")

    for num_steps in num_stepsList:
        start = time.time()
        tree = FinEquityBinomialTree()
        results = tree.value(
            stock_price,
            discount_curve,
            dividendCurve,
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

    payoff = FinEquityTreePayoffTypes.VANILLA_OPTION
    exercise = FinEquityTreeExerciseTypes.AMERICAN
    params = np.array([-1, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")

    for num_steps in num_stepsList:
        start = time.time()
        tree = FinEquityBinomialTree()
        results = tree.value(
            stock_price,
            discount_curve,
            dividendCurve,
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

    callOption = FinEquityVanillaOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
    delta = callOption.delta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    gamma = callOption.gamma(valuation_date, stock_price, discount_curve, dividendCurve, model)
    theta = callOption.theta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = FinEquityTreePayoffTypes.VANILLA_OPTION
    exercise = FinEquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")
    for num_steps in num_stepsList:
        start = time.time()
        tree = FinEquityBinomialTree()

        results = tree.value(
            stock_price,
            discount_curve,
            dividendCurve,
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

    payoff = FinEquityTreePayoffTypes.VANILLA_OPTION
    exercise = FinEquityTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strikePrice])

    testCases.header("NumSteps", "Results", "TIME")
    for num_steps in num_stepsList:
        start = time.time()
        tree = FinEquityBinomialTree()

        results = tree.value(
            stock_price,
            discount_curve,
            dividendCurve,
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
