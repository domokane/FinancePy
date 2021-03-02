###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityAmericanOption import FinEquityAmericanOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes, FinModelBlackScholesTypes
from financepy.utils.Date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def testFinEquityAmericanOption():

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(1, 1, 2017)
    stockPrice = 50.0
    interestRate = 0.06
    dividendYield = 0.04
    volatility = 0.40
    strikePrice = 50.0

    discount_curve = FinDiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = FinDiscountCurveFlat(valuation_date, dividendYield)

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = FinEquityAmericanOption(expiry_date, strikePrice, FinOptionTypes.EUROPEAN_PUT)

    model = FinModelBlackScholes(volatility, 
                                 FinModelBlackScholesTypes.CRR_TREE,
                                 100)

    value = putOption.value(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    delta = putOption.delta(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    gamma = putOption.gamma(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    theta = putOption.theta(valuation_date, stockPrice, discount_curve, dividendCurve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = FinEquityAmericanOption(expiry_date, strikePrice, FinOptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    numStepsList = [100, 200, 500, 1000, 2000]

    for numSteps in numStepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     numSteps)

        start = time.time()
        results = option.value(valuation_date, stockPrice, discount_curve, dividendCurve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", numSteps, results, duration)

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

    for numSteps in numStepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     numSteps)

        start = time.time()
        results = option.value(valuation_date, stockPrice, discount_curve, dividendCurve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", numSteps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = FinEquityAmericanOption(
        expiry_date,
        strikePrice,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    delta = callOption.delta(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    gamma = callOption.gamma(valuation_date, stockPrice, discount_curve, dividendCurve, model)
    theta = callOption.theta(valuation_date, stockPrice, discount_curve, dividendCurve, model)

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

    for numSteps in numStepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     numSteps)
        start = time.time()
        results = option.value(valuation_date, stockPrice, discount_curve, 
                               dividendCurve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_CALL_TREE", numSteps, results, duration)

    testCases.banner(
        "================== AMERICAN CALL =======================")
    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    option = FinEquityAmericanOption(expiry_date, strikePrice,
                                     FinOptionTypes.AMERICAN_CALL)

    for numSteps in numStepsList:

        model = FinModelBlackScholes(volatility,
                                     FinModelBlackScholesTypes.CRR_TREE,
                                     numSteps)

        start = time.time()

        results = option.value(valuation_date, stockPrice, discount_curve, 
                               dividendCurve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_CALL", numSteps, results, duration)

#    FinTest.TestReport(filename)

###############################################################################


testFinEquityAmericanOption()
testCases.compareTestCases()
