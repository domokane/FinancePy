# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


import add_fp_to_path

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_compound_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt1 = Date(1, 1, 2017)
    expiry_dt2 = Date(1, 1, 2018)
    k1 = 5.0
    k2 = 95.0
    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    num_steps_list = [100, 200, 500, 1000]

    stock_price = 85.0

    test_cases.header(
        "TYPE1", "TYPE2", "k_1", "k_2", "S", "TreeSteps", "Exact", "TreeValue"
    )

    for opt_type1 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        for opt_type2 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )

            for num_steps in num_steps_list:

                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                values_tree = cmpd_option.value_tree(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                test_cases.print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    num_steps,
                    value,
                    values_tree[0],
                )

    stock_price = 85.0

    test_cases.header(
        "TYPE1", "TYPE2", "k_1", "k_2", "S", "TreeSteps", "Exact", "TreeValue"
    )

    for opt_type1 in [OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT]:
        for opt_type2 in [OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )

            for num_steps in num_steps_list:

                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                values_tree = cmpd_option.value_tree(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                    num_steps,
                )

                test_cases.print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    num_steps,
                    value,
                    values_tree[0],
                )

    test_cases.header(
        "TYPE1",
        "TYPE2",
        "k_1",
        "k_2",
        "S",
        "Exact",
        "TreeSteps",
        "TreeValue",
        "Diff",
        "DELTA",
        "GAMMA",
        "THETA",
    )

    for opt_type1 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:
        for opt_type2 in [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]:

            cmpd_option = EquityCompoundOption(
                expiry_dt1, opt_type1, k1, expiry_dt2, opt_type2, k2
            )
            stock_prices = range(70, 100, 10)

            for stock_price in stock_prices:
                value = cmpd_option.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                delta = cmpd_option.delta(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                vega = cmpd_option.vega(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
                theta = cmpd_option.theta(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                values = cmpd_option.value_tree(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )

                diff = value - values[0]

                test_cases.print(
                    opt_type1,
                    opt_type2,
                    k1,
                    k2,
                    stock_price,
                    value,
                    num_steps,
                    values[0],
                    diff,
                    delta,
                    vega,
                    theta,
                )


########################################################################################

test_equity_compound_option()
test_cases.compare_test_cases()
