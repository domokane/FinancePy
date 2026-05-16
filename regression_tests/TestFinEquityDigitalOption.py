# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_digital_option import (
    EquityDigitalOption,
    FinDigitalOptionTypes,
)
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_digital_option():

    underlying_type = FinDigitalOptionTypes.CASH_OR_NOTHING

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    model = BlackScholes(volatility)
    import time

    call_option_values = []
    call_option_values_mc = []
    num_paths_list = [10000, 20000, 40000, 80000]

    #  [160000, 320000, 640000, 1280000, 2560000]

    test_cases.header("NumLoops", "ValueBS", "ValueMC", "TIME")

    for num_paths in num_paths_list:

        call_option = EquityDigitalOption(
            expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type
        )
        value = call_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        start = time.time()
        value_mc = call_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
        )
        end = time.time()
        duration = end - start
        test_cases.print(num_paths, value, value_mc, duration)

        call_option_values.append(value)
        call_option_values_mc.append(value_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(num_paths_list, call_option_values, color = 'b', label="Call Option")
    #    plt.plot(num_paths_list, call_option_values_mc, color = 'r', label = "Call Option MC")
    #    plt.xlabel("Num Loops")
    #    plt.legend(loc='best')

    stock_prices = range(50, 150, 50)
    call_option_values = []
    call_option_deltas = []
    call_option_vegas = []
    call_option_thetas = []

    for stock_price in stock_prices:
        call_option = EquityDigitalOption(
            expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type
        )
        value = call_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        delta = call_option.delta(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        vega = call_option.vega(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        theta = call_option.theta(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        call_option_values.append(value)
        call_option_deltas.append(delta)
        call_option_vegas.append(vega)
        call_option_thetas.append(theta)

    put_option_values = []
    put_option_deltas = []
    put_option_vegas = []
    put_option_thetas = []

    for stock_price in stock_prices:
        put_option = EquityDigitalOption(
            expiry_dt, 100.0, OptionTypes.EUROPEAN_PUT, underlying_type
        )
        value = put_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        delta = put_option.delta(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        vega = put_option.vega(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        theta = put_option.theta(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        put_option_values.append(value)
        put_option_deltas.append(delta)
        put_option_vegas.append(vega)
        put_option_thetas.append(theta)


########################################################################################

test_equity_digital_option()
test_cases.compare_test_cases()
