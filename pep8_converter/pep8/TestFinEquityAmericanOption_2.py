# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
from financepy.products.equity.equity_american_option import (
    EquityAmericanOption,
)
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)



########################################################################################


def test_equity_american_option():

    value_dt = Date(1, 1, 2016)
    expiry_dt = Date(1, 1, 2017)
    stock_price = 50.0
    interest_rate = 0.06
    dividend_yield = 0.04
    volatility = 0.40
    strike_price = 50.0

    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    test_cases.banner("================== EUROPEAN PUT =======================")

    put_option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT)

    num_steps = 4

    model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)

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

    test_cases.header("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
    test_cases.print("EUROPEAN_PUT_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT)

    test_cases.header("opt_type", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    num_steps_list = [100, 200, 1000]

    for num_steps in num_steps_list:

        model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)

        start = time.time()
        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        end = time.time()
        duration = end - start
        test_cases.print("EUROPEAN_PUT_TREE", num_steps, results, duration)

    test_cases.banner("================== AMERICAN PUT =======================")

    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.AMERICAN_PUT)

    test_cases.header("opt_type", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)

        start = time.time()
        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        end = time.time()
        duration = end - start
        test_cases.print("AMERICAN_PUT", num_steps, results, duration)

    test_cases.banner("================== EUROPEAN CALL =======================")

    call_option = EquityAmericanOption(
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

    test_cases.header("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
    test_cases.print("EUROPEAN_CALL_BS", value, delta, gamma, theta)

    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL)

    test_cases.header("opt_type", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    for num_steps in num_steps_list:

        model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)
        start = time.time()
        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        end = time.time()
        duration = end - start
        test_cases.print("EUROPEAN_CALL_TREE", num_steps, results, duration)

    test_cases.banner("================== AMERICAN CALL =======================")
    test_cases.header("opt_type", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.AMERICAN_CALL)

    for num_steps in num_steps_list:

        model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)

        start = time.time()

        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        end = time.time()
        duration = end - start
        test_cases.print("AMERICAN_TREE_CALL", num_steps, results, duration)

    test_cases.banner("================== AMERICAN PUT =======================")
    test_cases.header("opt_type", "NUMSTEPS", "VALUE DELTA GAMMA THETA", "TIME")

    option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.AMERICAN_PUT)

    for num_steps in num_steps_list:

        model = BlackScholes(volatility, BlackScholesTypes.CRR_TREE, num_steps)

        start = time.time()

        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        end = time.time()
        duration = end - start
        test_cases.print("AMERICAN_TREE_PUT", num_steps, results, duration)


#    fin_test.test_report(filename)

########################################################################################


test_equity_american_option()
test_cases.compare_test_cases()

########################################################################################

