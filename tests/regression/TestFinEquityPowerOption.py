# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 14:03:26 2026

@author: Dominic
"""

# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time

import add_fp_to_path

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_power_option import EquityPowerOption
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)


########################################################################################


def test_equity_power_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 7, 2015)

    stock_price = 100.0
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    power = 2

    # For power=2, K=10000 corresponds to an exercise boundary
    # S(T) = sqrt(K) = 100.
    strike_price = 10000.0

    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    ###########################################################################
    # Convergence of MC to analytical value
    ###########################################################################

    num_paths_list = [10000, 20000, 40000]

    test_cases.header(
        "NUMPATHS",
        "VALUE_BS",
        "VALUE_MC",
        "TIME",
    )

    call_option = EquityPowerOption(
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
        power,
    )

    for num_paths in num_paths_list:

        value = call_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        start = time.time()

        value_mc = call_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
        ).value

        end = time.time()

        test_cases.print(
            num_paths,
            value,
            value_mc,
            end - start,
        )

    ###########################################################################
    # Call values across stock prices
    ###########################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000
    seed = 1967

    test_cases.header(
        "STOCK_PRICE",
        "NUMPATHS",
        "CALL_VALUE_BS",
        "CALL_VALUE_MC",
        "TIME",
    )

    for stock_price in stock_prices:

        call_option = EquityPowerOption(
            expiry_dt,
            strike_price,
            OptionTypes.EUROPEAN_CALL,
            power,
        )

        value = call_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        start = time.time()

        value_mc = call_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
        ).value

        end = time.time()

        test_cases.print(
            stock_price,
            num_paths,
            value,
            value_mc,
            end - start,
        )

    ###########################################################################
    # Put values across stock prices
    ###########################################################################

    test_cases.header(
        "STOCK_PRICE",
        "NUMPATHS",
        "PUT_VALUE_BS",
        "PUT_VALUE_MC",
        "TIME",
    )

    for stock_price in stock_prices:

        put_option = EquityPowerOption(
            expiry_dt,
            strike_price,
            OptionTypes.EUROPEAN_PUT,
            power,
        )

        value = put_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        start = time.time()

        value_mc = put_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
        ).value

        end = time.time()

        test_cases.print(
            stock_price,
            num_paths,
            value,
            value_mc,
            end - start,
        )

    ###########################################################################
    # n = 1 should reproduce vanilla Black-Scholes
    ###########################################################################

    test_cases.header(
        "STOCK_PRICE",
        "POWER_CALL",
        "VANILLA_CALL",
        "DIFFERENCE",
    )

    strike = 100.0

    for stock_price in stock_prices:

        power_option = EquityPowerOption(
            expiry_dt,
            strike,
            OptionTypes.EUROPEAN_CALL,
            1,
        )

        vanilla_option = EquityVanillaOption(
            expiry_dt,
            strike,
            OptionTypes.EUROPEAN_CALL,
        )

        power_value = power_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        vanilla_value = vanilla_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        test_cases.print(
            stock_price,
            power_value,
            vanilla_value,
            power_value - vanilla_value,
        )


########################################################################################


test_equity_power_option()
test_cases.compare_test_cases()
