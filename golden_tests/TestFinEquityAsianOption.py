# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import matplotlib.pyplot as plt

import glob
from os.path import dirname, basename, join
import traceback

import add_fp_to_path

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.equity.equity_asian_option import AsianOptionValuationMethods
from financepy.products.equity.equity_asian_option import EquityAsianOption
from financepy.utils.global_types import OptionTypes

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

PLOT_FLAG = False

test_convergence = False
test_time_evolution = False
test_mc_timings = True

########################################################################################


def test_convergence_fn():

    value_dt = Date(1, 1, 2014)
    start_averaging_date = Date(1, 6, 2014)
    expiry_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 120  # daily as we have a half year
    accrued_avg = None
    k = 100
    seed = 1976

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    asian_option = EquityAsianOption(
        start_averaging_date,
        expiry_dt,
        k,
        OptionTypes.EUROPEAN_CALL,
        num_observations,
    )

    test_cases.header(
        "K", "Geometric", "Turnbull_Wakeman", "Curran", "FastMC", "FastMC_CV"
    )

    values_turnbull = []
    values_curran = []
    values_geometric = []
    values_mc_fast = []
    values_mc_cv = []

    num_paths_list = [5000]

    for num_paths in num_paths_list:

        accrued_avg = stock_price * 1.1

        value_mc_fast = asian_option.value_mc_fast(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        value_mc_cv = asian_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        value_geometric = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.GEOMETRIC,
            accrued_avg,
        )

        value_turnbull_wakeman = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.TURNBULL_WAKEMAN,
            accrued_avg,
        )

        value_curran = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.CURRAN,
            accrued_avg,
        )

        values_geometric.append(value_geometric)
        values_turnbull.append(value_turnbull_wakeman)
        values_curran.append(value_curran)
        values_mc_fast.append(value_mc_fast)
        values_mc_cv.append(value_mc_cv)

        test_cases.print(
            num_paths,
            value_geometric,
            value_turnbull_wakeman,
            value_curran,
            value_mc_fast,
            value_mc_cv,
        )


if PLOT_FLAG:
    import matplotlib.pyplot as plt

    x = num_paths_list
    plt.figure(figsize=(8, 6))
    plt.plot(x, values_geometric, label="Geometric")
    plt.plot(x, values_turnbull, label="Turbull_Wakeman")
    plt.plot(x, values_curran, label="Curran")
    plt.plot(x, values_mc_fast, label="MC_Fast")
    plt.plot(x, values_mc_cv, label="MC_CV")
    plt.legend()
    plt.xlabel("Number of Paths")
    plt.show()

########################################################################################


def test_time_evolution_fn():

    start_averaging_date = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 100  # weekly as we have a year
    accrued_avg = None
    k = 100
    seed = 1976

    model = BlackScholes(volatility)

    asian_option = EquityAsianOption(
        start_averaging_date,
        expiry_dt,
        k,
        OptionTypes.EUROPEAN_CALL,
        num_observations,
    )

    test_cases.header(
        "Date",
        "Geometric",
        "Turnbull_Wakeman",
        "Curran",
        "FastMC",
        "FastMC_CV",
    )

    values_turnbull = []
    values_curran = []
    values_geometric = []
    values_mc_fast = []
    values_mc_cv = []

    value_dts = []
    value_dts.append(Date(1, 4, 2014))
    value_dts.append(Date(1, 6, 2014))
    value_dts.append(Date(1, 8, 2014))
    value_dts.append(Date(1, 2, 2015))
    value_dts.append(Date(1, 4, 2015))
    value_dts.append(Date(1, 6, 2015))
    value_dts.append(Date(1, 8, 2015))

    num_paths = 10000

    for value_dt in value_dts:

        accrued_avg = stock_price * 0.9

        discount_curve = DiscountCurveFlat(value_dt, interest_rate)
        dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

        value_mc_fast = asian_option.value_mc_fast(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        value_mc_cv = asian_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        value_geometric = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.GEOMETRIC,
            accrued_avg,
        )

        value_turnbull_wakeman = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.TURNBULL_WAKEMAN,
            accrued_avg,
        )

        value_curran = asian_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            AsianOptionValuationMethods.CURRAN,
            accrued_avg,
        )

        values_geometric.append(value_geometric)
        values_turnbull.append(value_turnbull_wakeman)
        values_curran.append(value_curran)
        values_mc_fast.append(value_mc_fast)
        values_mc_cv.append(value_mc_cv)

        test_cases.print(
            str(value_dt),
            value_geometric,
            value_turnbull_wakeman,
            value_curran,
            value_mc_fast,
            value_mc_cv,
        )


if PLOT_FLAG:

    x = [dt.date() for dt in value_dts]
    plt.figure(figsize=(8, 6))
    plt.plot(x, values_geometric, label="Geometric")
    plt.plot(x, values_turnbull, label="Turbull_Wakeman")
    plt.plot(x, values_curran, label="Curran")
    plt.plot(x, values_mc_fast, label="MC_Fast")
    plt.plot(x, values_mc_cv, label="MC_CV")
    plt.legend()
    plt.xlabel("Valuation Date")
    plt.show()

########################################################################################


def test_mc_timings_fn():

    value_dt = Date(1, 1, 2014)
    start_averaging_date = Date(1, 6, 2014)
    expiry_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.30
    dividend_yield = 0.10
    num_observations = 120  # daily as we have a half year
    accrued_avg = None
    k = 100
    seed = 1976

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    asian_option = EquityAsianOption(
        start_averaging_date,
        expiry_dt,
        k,
        OptionTypes.EUROPEAN_CALL,
        num_observations,
    )

    test_cases.header(
        "NUMPATHS", "VALUE", "TIME", "VALUE_MC", "TIME", "VALUE_MC_CV", "TIME"
    )

    values_mc = []
    values_mc_fast = []
    values_mc_fast_cv = []

    tvalues_mc = []
    tvalues_mc_fast = []
    tvalues_mc_fast_cv = []

    num_paths_list = [5000]

    for num_paths in num_paths_list:

        accrued_avg = stock_price * 1.1

        start = time.time()
        value_mc = asian_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        end = time.time()
        t_mc = end - start

        start = time.time()
        value_mc_fast = asian_option.value_mc_fast(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        end = time.time()
        t_mc_fast = end - start

        start = time.time()
        value_mc_fast_cv = asian_option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths,
            seed,
            accrued_avg,
        )

        end = time.time()
        t_mc_fast_cv = end - start

        values_mc.append(value_mc)
        values_mc_fast.append(value_mc_fast)
        values_mc_fast_cv.append(value_mc_fast_cv)

        tvalues_mc.append(t_mc)
        tvalues_mc_fast.append(t_mc_fast)
        tvalues_mc_fast_cv.append(t_mc_fast_cv)

        test_cases.print(
            num_paths,
            value_mc,
            t_mc,
            value_mc_fast,
            t_mc_fast,
            value_mc_fast_cv,
            t_mc_fast_cv,
        )


########################################################################################

if PLOT_FLAG:

    x = num_paths_list
    plt.figure(figsize=(8, 6))
    plt.plot(x, values_mc, label="Basic MC")
    plt.plot(x, values_mc_fast, label="MC_Fast")
    plt.plot(x, values_mc_fast_cv, label="MC_Fast CV")
    plt.legend()
    plt.xlabel("Number of Paths")
    plt.show()


test_convergence_fn()
test_mc_timings_fn()
test_time_evolution_fn()

test_cases.compare_test_cases()
