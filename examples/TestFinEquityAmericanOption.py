###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityAmericanOption import FinEquityAmericanOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.market.curves.FinDiscountCurveFlat import DiscountCurveFlat
from financepy.models.black_scholes import FinModelBlackScholes, FinModelBlackScholesTypes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def testFinEquityAmericanOption():

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    stock_price = 50.0
    interestRate = 0.06
    dividendYield = 0.04
    volatility = 0.40
    strikePrice = 50.0

    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = FinEquityAmericanOption(expiry_date, strikePrice, FinOptionTypes.EUROPEAN_PUT)

    model = FinModelBlackScholes(volatility, 
                                 FinModelBlackScholesTypes.CRR_TREE,
                                 100)

    value = putOption.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
    delta = putOption.delta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    gamma = putOption.gamma(valuation_date, stock_price, discount_curve, dividendCurve, model)
    theta = putOption.theta(valuation_date, stock_price, discount_curve, dividendCurve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = FinEquityAmericanOption(expiry_date, strikePrice, FinOptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    num_stepsList = [100, 200, 500, 1000, 2000]

    for num_steps in num_stepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", num_steps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = FinEquityAmericanOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_stepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", num_steps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = FinEquityAmericanOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valuation_date, stock_price, discount_curve, dividendCurve, model)
    delta = callOption.delta(valuation_date, stock_price, discount_curve, dividendCurve, model)
    gamma = callOption.gamma(valuation_date, stock_price, discount_curve, dividendCurve, model)
    theta = callOption.theta(valuation_date, stock_price, discount_curve, dividendCurve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = FinEquityAmericanOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_stepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     num_steps)
        start = time.time()
        results = option.value(valuation_date, stock_price, discount_curve,
                               dividendCurve, model)
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

    option = FinEquityAmericanOption(expiry_date, strikePrice,
                                     FinOptionTypes.AMERICAN_CALL)

    for num_steps in num_stepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     num_steps)

        start = time.time()

        results = option.value(valuation_date, stock_price, discount_curve,
                               dividendCurve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_CALL", num_steps, results, duration)

#    FinTest.TestReport(filename)

###############################################################################


testFinEquityAmericanOption()
testCases.compareTestCases()
