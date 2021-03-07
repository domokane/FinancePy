###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import FinOptionTypes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes, FinModelBlackScholesTypes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def testEquityAmericanOption():

    valueDate = FinDate(1, 1, 2016)
    expiryDate = FinDate(1, 1, 2017)
    stockPrice = 50.0
    interestRate = 0.06
    dividendYield = 0.04
    volatility = 0.40
    strike_price = 50.0

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    testCases.banner("================== EUROPEAN PUT =======================")

    putOption = EquityAmericanOption(expiryDate, strike_price, FinOptionTypes.EUROPEAN_PUT)

    model = BlackScholes(volatility,
                         FinModelBlackScholesTypes.CRR_TREE,
                         100)

    value = putOption.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
    delta = putOption.delta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    gamma = putOption.gamma(valueDate, stockPrice, discountCurve, dividend_curve, model)
    theta = putOption.theta(valueDate, stockPrice, discountCurve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(expiryDate, strike_price, FinOptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    numStepsList = [100, 200, 500, 1000, 2000]

    for numSteps in numStepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             numSteps)

        start = time.time()
        results = option.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", numSteps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = EquityAmericanOption(
        expiryDate,
        strike_price,
        FinOptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for numSteps in numStepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             numSteps)

        start = time.time()
        results = option.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", numSteps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    callOption = EquityAmericanOption(
        expiryDate,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)
    value = callOption.value(valueDate, stockPrice, discountCurve, dividend_curve, model)
    delta = callOption.delta(valueDate, stockPrice, discountCurve, dividend_curve, model)
    gamma = callOption.gamma(valueDate, stockPrice, discountCurve, dividend_curve, model)
    theta = callOption.theta(valueDate, stockPrice, discountCurve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiryDate,
        strike_price,
        FinOptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for numSteps in numStepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             numSteps)
        start = time.time()
        results = option.value(valueDate, stockPrice, discountCurve, 
                               dividend_curve, model)
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

    option = EquityAmericanOption(expiryDate, strike_price,
                                     FinOptionTypes.AMERICAN_CALL)

    for numSteps in numStepsList:

        model = BlackScholes(volatility,
                             FinModelBlackScholesTypes.CRR_TREE,
                             numSteps)

        start = time.time()

        results = option.value(valueDate, stockPrice, discountCurve, 
                               dividend_curve, model)

        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_CALL", numSteps, results, duration)

#    FinTest.TestReport(filename)

###############################################################################


testEquityAmericanOption()
testCases.compareTestCases()
