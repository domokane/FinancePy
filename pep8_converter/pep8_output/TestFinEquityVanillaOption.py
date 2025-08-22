# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
import numpy as np
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.error import FinError
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test__equity_vanilla_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    num_paths_list = [10000, 20000, 40000]  # , 80000, 160000, 320000]

    test_cases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for num_paths in num_paths_list:

        call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)

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


    stock_prices = range(80, 120, 10)
    num_paths = 100000

    test_cases.header(
        "NUMPATHS",
        "CALL_VALUE_BS",
        "CALL_VALUE_MC",
        "CALL_VALUE_MC_SOBOL",
        "TIME",
    )
    use_sobol = True

    for stock_price in stock_prices:

        call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)

        value = call_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        start = time.time()

        use_sobol = False
        value_mc1 = call_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            use_sobol,
        )

        use_sobol = True
        value_mc2 = call_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            use_sobol,
        )

        end = time.time()
        duration = end - start
        test_cases.print(num_paths, value, value_mc1, value_mc2, duration)


    stock_prices = range(80, 120, 10)
    num_paths = 100000

    test_cases.header(
        "NUMPATHS",
        "PUT_VALUE_BS",
        "PUT_VALUE_MC",
        "PUT_VALUE_MC_SOBOL",
        "TIME",
    )

    for stock_price in stock_prices:

        put_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_PUT)

        value = put_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        start = time.time()

        use_sobol = False
        value_mc1 = put_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            use_sobol,
        )

        use_sobol = True
        value_mc2 = put_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            use_sobol,
        )

        end = time.time()
        duration = end - start
        test_cases.print(num_paths, value, value_mc1, value_mc2, duration)


    stock_prices = range(80, 120, 10)

    test_cases.header(
        "STOCK PRICE",
        "CALL_VALUE_BS",
        "CALL_DELTA_BS",
        "CALL_VEGA_BS",
        "CALL_THETA_BS",
        "CALL_RHO_BS",
        "CALL_VANNA_BS",
    )

    for stock_price in stock_prices:

        call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)
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
        rho = call_option.rho(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        vanna = call_option.vanna(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        test_cases.print(stock_price, value, delta, vega, theta, rho, vanna)


    test_cases.header(
        "STOCK PRICE",
        "PUT_VALUE_BS",
        "PUT_DELTA_BS",
        "PUT_VEGA_BS",
        "PUT_THETA_BS",
        "PUT_RHO_BS",
        "PUT_VANNA_BS",
    )

    for stock_price in stock_prices:

        put_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_PUT)

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
        rho = put_option.rho(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        vanna = put_option.vanna(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )
        test_cases.print(stock_price, value, delta, vega, theta, rho, vanna)

################################################################################


def test_implied_volatility_new():

    value_dt = Date(1, 1, 2015)
    stock_price = 100.0
    interest_rate = 0.05
    dividend_yield = 0.03
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    strikes = np.linspace(50, 150, 11)
    times_to_expiry = [0.003, 0.01, 0.1, 0.5, 1.0]
    sigmas = np.arange(1, 100, 5) / 100.0
    opt_types = [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]

    test_cases.header(
        "OPT_TYPE",
        "TEXP",
        "STOCK_PRICE",
        "STRIKE",
        "INTRINSIC",
        "VALUE",
        "INPUT_VOL",
        "IMPLIED_VOL",
    )

    tol = 1e-5
    num_tests = 0
    num_fails = 0

    for vol in sigmas:

        model = BlackScholes(vol)

        for time_to_expiry in times_to_expiry:

            expiry_dt = value_dt.add_years(time_to_expiry)

            for strike in strikes:

                for opt_type in opt_types:

                    option = EquityVanillaOption(expiry_dt, strike, opt_type)

                    value = option.value(
                        value_dt,
                        stock_price,
                        discount_curve,
                        dividend_curve,
                        model,
                    )

                    intrinsic = option.intrinsic(
                        value_dt, stock_price, discount_curve, dividend_curve
                    )

                    # I remove the cases where the time value is zero
                    # This is arbitrary but 1e-10 seems good enough to me

                    implied_vol = -999

                    if value - intrinsic > 1e-10:

                        implied_vol = option.implied_volatility(
                            value_dt,
                            stock_price,
                            discount_curve,
                            dividend_curve,
                            value,
                        )

                    num_tests += 1

                    err_vol = np.abs(implied_vol - vol)

                    if err_vol > tol:

                        test_cases.print(
                            opt_type,
                            time_to_expiry,
                            stock_price,
                            strike,
                            intrinsic,
                            value,
                            vol,
                            implied_vol,
                        )

                        # These fails include ones due to the zero time value
                        num_fails += 1

                        test_cases.print(
                            opt_type,
                            time_to_expiry,
                            stock_price,
                            strike,
                            stock_price,
                            value,
                            vol,
                            implied_vol,
                        )


################################################################################

#    print("Num Tests", num_tests, "num_fails", num_fails)



if 1 == 0:
    value_dt = Date(30, 11, 2021)
    expiry_dt = value_dt.add_years(1)

    stock_price = 100
    volatility = 0.20
    model = BlackScholes(volatility)

    discount_curve = DiscountCurveFlat(value_dt, 0.05)
    dividend_curve = DiscountCurveFlat(value_dt, 0.0)

    call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)

    value = call_option.value(value_dt, 105.0, discount_curve, dividend_curve, model)

################################################################################

else:
    test__equity_vanilla_option()

    start = time.time()
    test_implied_volatility_new()
    end = time.time()
    elapsed = end - start

    # print("Elapsed:", elapsed)

    test_cases.compare_test_cases()
