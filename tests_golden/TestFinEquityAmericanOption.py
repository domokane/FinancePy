###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import FinOptionTypes
from financepy.market.curves.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes, FinModelBlackScholesTypes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def testEquityAmericanOption():

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    stock_price = 50.0
    interest_rate = 0.06
    dividend_yield = 0.04
    volatility = 0.40
    strike_price = 50.0

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = EquityAmericanOption(expiry_date, strike_price, FinOptionTypes.EUROPEAN_PUT)

    model = BlackScholes(volatility,
                         FinModelBlackScholesTypes.CRR_TREE,
                         100)

    value = putOption.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
    delta = putOption.delta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    gamma = putOption.gamma(valuation_date, stock_price, discount_curve, dividend_curve, model)
    theta = putOption.theta(valuation_date, stock_price, discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(expiry_date, strike_price, FinOptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    num_stepsList = [100, 200, 500, 1000, 2000]

    for num_steps in num_stepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", num_steps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        FinOptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_stepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", num_steps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = EquityAmericanOption(
        expiry_date,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valuation_date, stock_price, discount_curve, dividend_curve, model)
    delta = callOption.delta(valuation_date, stock_price, discount_curve, dividend_curve, model)
    gamma = callOption.gamma(valuation_date, stock_price, discount_curve, dividend_curve, model)
    theta = callOption.theta(valuation_date, stock_price, discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_stepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             num_steps)
        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_CALL_TREE", num_steps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")
    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    option = EquityAmericanOption(expiry_date, strike_price,
                                     FinOptionTypes.AMERICAN_CALL)

    for num_steps in num_stepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_CALL", num_steps, results, duration)

#    FinTest.TestReport(filename)

###############################################################################


testEquityAmericanOption()
testCases.compareTestCases()
