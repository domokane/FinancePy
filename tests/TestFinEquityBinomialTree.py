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

    stockPrice = 50.0
    riskFreeRate = 0.06
    dividendYield = 0.04
    volatility = 0.40

    valueDate = FinDate(1, 1, 2016)
    expiryDate = FinDate(1, 1, 2017)

    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, riskFreeRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    numStepsList = [100, 500, 1000, 2000, 5000]

    strike_price = 50.0

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = EquityVanillaOption(
        expiryDate,
        strike_price,
        FinOptionTypes.EUROPEAN_PUT)
    value = putOption.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
    delta = putOption.delta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    gamma = putOption.gamma(valueDate, stockPrice, discountCurve, dividend_curve, model)
    theta = putOption.theta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([-1, strike_price])

    testCases.header("NumSteps", "Results", "TIME")

    for numSteps in numStepsList:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
            stockPrice,
            discountCurve,
            dividend_curve,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)
        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([-1, strike_price])

    testCases.header("NumSteps", "Results", "TIME")

    for numSteps in numStepsList:
        start = time.time()
        tree = EquityBinomialTree()
        results = tree.value(
            stockPrice,
            discountCurve,
            dividend_curve,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)
        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = EquityVanillaOption(
        expiryDate,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
    delta = callOption.delta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    gamma = callOption.gamma(valueDate, stockPrice, discountCurve, dividend_curve, model)
    theta = callOption.theta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    testCases.header("BS Value", "BS Delta", "BS Gamma", "BS Theta")
    testCases.print(value, delta, gamma, theta)

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.EUROPEAN
    params = np.array([1.0, strike_price])

    testCases.header("NumSteps", "Results", "TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
            stockPrice,
            discountCurve,
            dividend_curve,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)

        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    exercise = EquityTreeExerciseTypes.AMERICAN
    params = np.array([1.0, strike_price])

    testCases.header("NumSteps", "Results", "TIME")
    for numSteps in numStepsList:
        start = time.time()
        tree = EquityBinomialTree()

        results = tree.value(
            stockPrice,
            discountCurve,
            dividend_curve,
            volatility,
            numSteps,
            valueDate,
            payoff,
            expiryDate,
            payoff,
            exercise,
            params)

        end = time.time()
        duration = end - start
        testCases.print(numSteps, results, duration)

###############################################################################

test_FinBinomialTree()
testCases.compareTestCases()
