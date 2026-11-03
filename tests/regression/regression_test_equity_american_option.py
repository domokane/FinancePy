# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time

import add_fp_to_path

from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.models.bs_american_analytic import BlackScholesAlgorithmTypes

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

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    test_cases.banner("================== EUROPEAN PUT =======================")

    put_option = EquityAmericanOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT)

    num_steps = 4

    model = BlackScholes(volatility)

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

        model = BlackScholes(volatility)

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

        model = BlackScholes(volatility)

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

        model = BlackScholes(volatility)
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

        model = BlackScholes(volatility)

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

        model = BlackScholes(volatility)

        start = time.time()

        results = option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        end = time.time()
        duration = end - start
        test_cases.print("AMERICAN_TREE_PUT", num_steps, results, duration)


########################################################################################

def test_equity_american_option_algorithms():

    value_dt = Date(1, 1, 2016)
    expiry_dt = Date(1, 1, 2017)

    stock_price = 50.0
    strike_price = 50.0

    interest_rate = 0.06
    dividend_yield = 0.04
    volatility = 0.40

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    model = BlackScholes(volatility)

    ############################################################################
    # Algorithms
    ############################################################################

    european_algorithms = [
        BlackScholesAlgorithmTypes.DEFAULT,
        BlackScholesAlgorithmTypes.ANALYTICAL,
        BlackScholesAlgorithmTypes.CRR_TREE,
        BlackScholesAlgorithmTypes.LSMC,
        BlackScholesAlgorithmTypes.FINITE_DIFFERENCE,
    ]

    american_algorithms = [
        BlackScholesAlgorithmTypes.DEFAULT,
        BlackScholesAlgorithmTypes.CRR_TREE,
        BlackScholesAlgorithmTypes.BARONE_ADESI_WHALEY,
        BlackScholesAlgorithmTypes.LSMC,
        BlackScholesAlgorithmTypes.BJERKSUND_STENSLAND_1993,
        BlackScholesAlgorithmTypes.BJERKSUND_STENSLAND_2002,
        BlackScholesAlgorithmTypes.FINITE_DIFFERENCE,
    ]

    ############################################################################
    # Option types
    ############################################################################

    option_tests = [
        (
            "EUROPEAN_CALL",
            OptionTypes.EUROPEAN_CALL,
            european_algorithms,
        ),
        (
            "EUROPEAN_PUT",
            OptionTypes.EUROPEAN_PUT,
            european_algorithms,
        ),
        (
            "AMERICAN_CALL",
            OptionTypes.AMERICAN_CALL,
            american_algorithms,
        ),
        (
            "AMERICAN_PUT",
            OptionTypes.AMERICAN_PUT,
            american_algorithms,
        ),
    ]

    ############################################################################
    # Run tests
    ############################################################################

    for option_name, option_type, algorithms in option_tests:

        test_cases.banner(
            f"================== {option_name} ======================="
        )

        test_cases.header(
            "OPTION",
            "ALGORITHM",
            "VALUE",
            "TIME",
        )

        option = EquityAmericanOption(
            expiry_dt,
            strike_price,
            option_type,
        )

        for algo_type in algorithms:

            start = time.time()

            value = option.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                algo_type=algo_type,
            )

            end = time.time()

            test_cases.print(
                option_name,
                algo_type.name,
                value,
                end - start,
            )


########################################################################################

#    FinTest.TestReport(filename)


test_equity_american_option()
test_equity_american_option_algorithms()
test_cases.compare_test_cases()
