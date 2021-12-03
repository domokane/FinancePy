###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


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

    put_option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_PUT)

    model = BlackScholes(volatility,
                         BlackScholesTypes.CRR_TREE,
                         100)

    value = put_option.value(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    delta = put_option.delta(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    gamma = put_option.gamma(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)
    theta = put_option.theta(valuation_date, stock_price,
                             discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_PUT)

    testCases.header("OPTION_TYPE", "NUMSTEPS",
                     "VALUE DELTA GAMMA THETA", "TIME")

    num_steps_list = [100, 200, 500, 1000, 2000]

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price,
                               discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("EUROPEAN_PUT_TREE", num_steps, results, duration)

    testCases.banner("================== AMERICAN PUT =======================")

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.AMERICAN_PUT)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
                             num_steps)

        start = time.time()
        results = option.value(valuation_date, stock_price,
                               discount_curve, dividend_curve, model)
        end = time.time()
        duration = end - start
        testCases.print("AMERICAN_PUT", num_steps, results, duration)

    testCases.banner(
        "================== EUROPEAN CALL =======================")

    call_option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.EUROPEAN_CALL)

    value = call_option.value(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    delta = call_option.delta(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    gamma = call_option.gamma(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)
    theta = call_option.theta(valuation_date, stock_price,
                              discount_curve, dividend_curve, model)

    testCases.header("OPTION_TYPE", "VALUE", "DELTA", "GAMMA", "THETA")
    testCases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(
        expiry_date,
        strike_price,
        OptionTypes.EUROPEAN_CALL)

    testCases.header(
        "OPTION_TYPE",
        "NUMSTEPS",
        "VALUE DELTA GAMMA THETA",
        "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
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
                                  OptionTypes.AMERICAN_CALL)

    for num_steps in num_steps_list:

        model = BlackScholes(volatility,
                             BlackScholesTypes.CRR_TREE,
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
