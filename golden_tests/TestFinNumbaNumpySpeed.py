# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import matplotlib.pyplot as plt
import numpy as np

import add_fp_to_path

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_numba_numpy_speed(use_sobol):

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 1999

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)

    use_sobol_int = int(use_sobol)

    test_cases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)

    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_yield, model
    )

    num_points = 20
    v_exact = [value] * num_points

    # DO UP TO 100K AS IT IS SLOW

    num_paths_list = np.arange(1, num_points + 1, 1) * 100000

    nonumba_nonumpy_v = []
    nonumba_nonumpy_t = []

    print("PURE PYTHON")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_nonumba_nonumpy(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        nonumba_nonumpy_v.append(value_mc)
        nonumba_nonumpy_t.append(duration + 1e-10)

    numpy_only_v = []
    numpy_only_t = []

    print("NUMPY ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_only(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numpy_only_v.append(value_mc)
        numpy_only_t.append(duration + 1e-10)

    #    speedUp = np.array(nonumba_nonumpy_t)/np.array(numpy_only_t)
    #    print(numpy_only_t)
    #    print(nonumba_nonumpy_t)
    #    print(speedUp)

    if use_sobol:
        title = "SOBOL: PURE PYTHON VS NUMPY"
    else:
        title = "PSEUDORANDOM: PURE PYTHON VS NUMPY"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, nonumba_nonumpy_t, "o-", label="PURE PYTHON")
    plt.plot(num_paths_list, numpy_only_t, "o-", label="NUMPY ONLY")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, nonumba_nonumpy_v, "o-", label="PURE PYTHON")
    plt.plot(num_paths_list, numpy_only_v, "o-", label="NUMPY ONLY")

    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)

    # DO UP TO 10 MILLION NOW THAT WE HAVE NUMPY

    num_paths_list = np.arange(1, num_points + 1, 1) * 1000000

    numpy_only_v = []
    numpy_only_t = []

    print("NUMPY ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_only(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numpy_only_v.append(value_mc)
        numpy_only_t.append(duration)

    numba_numpy_v = []
    numba_numpy_t = []

    print("NUMBA+NUMPY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_numba(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numba_numpy_v.append(value_mc)
        numba_numpy_t.append(duration)

    numba_only_v = []
    numba_only_t = []

    print("NUMBA ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_only(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numba_only_v.append(value_mc)
        numba_only_t.append(duration)

    numba_parallel_v = []
    numba_parallel_t = []

    print("NUMBA PARALLEL")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_parallel(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numba_parallel_v.append(value_mc)
        numba_parallel_t.append(duration)

    #    speedUp = np.array(numba_only_t)/np.array(numba_parallel_t)
    #    print("PARALLEL:", speedUp)

    # COMPUTED USING NUMSTEPS FROM 1M to 10M

    cpp_t = np.array(
        [
            0.075,
            0.155,
            0.223,
            0.313,
            0.359,
            0.421,
            0.495,
            0.556,
            0.64,
            0.702,
            0.765,
            0.841,
            0.923,
            0.982,
            1.05,
            1.125,
            1.195,
            1.261,
            1.333,
            1.408,
        ]
    )

    cpp_v = np.array(
        [
            9.30872,
            9.29576,
            9.29422,
            9.29832,
            9.29863,
            9.30153,
            9.2994,
            9.3025,
            9.29653,
            9.29875,
            9.29897,
            9.29996,
            9.29931,
            9.29796,
            9.29784,
            9.2992,
            9.3001,
            9.30093,
            9.29876,
            9.29921,
        ]
    )

    if use_sobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, numpy_only_t, "o-", label="NUMPY ONLY")
    plt.plot(num_paths_list, numba_numpy_t, "o-", label="NUMBA + NUMPY")

    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    if use_sobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, numpy_only_t, "o-", label="NUMPY ONLY")
    plt.plot(num_paths_list, numba_numpy_t, "o-", label="NUMBA + NUMPY")
    plt.plot(num_paths_list, numba_only_t, "o-", label="NUMBA ONLY")
    plt.plot(num_paths_list, numba_parallel_t, "o-", label="NUMBA PARALLEL")

    if use_sobol is False:
        plt.plot(num_paths_list, cpp_t, "o-", label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, numba_only_v, "o-", label="NUMBA ONLY")
    plt.plot(num_paths_list, cpp_v, "o-", label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)


########################################################################################


def test_fin_numba_numba_parallel(use_sobol):

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 2021

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)

    use_sobol_int = int(use_sobol)

    test_cases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    call_option = EquityVanillaOption(expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL)

    value = call_option.value(
        value_dt, stock_price, discount_curve, dividend_yield, model
    )

    num_points = 20
    v_exact = [value] * num_points

    num_paths_list = np.arange(1, num_points + 1, 1) * 1000000

    numba_only_v = []
    numba_only_t = []

    print("NUMBA ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_only(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numba_only_v.append(value_mc)
        numba_only_t.append(duration)

    numba_parallel_v = []
    numba_parallel_t = []

    print("NUMBA PARALLEL")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_parallel(
            value_dt,
            stock_price,
            discount_curve,
            dividend_yield,
            model,
            num_paths,
            seed,
            use_sobol_int,
        )
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (num_paths, value, value_mc, duration))

        numba_parallel_v.append(value_mc)
        numba_parallel_t.append(duration)

    import matplotlib.pyplot as plt

    if use_sobol:
        title = "SOBOL: NUMBA VS NUMBA + PARALLEL"
    else:
        title = "PSEUDORANDOM: NUMBA VS NUMBA + PARALLEL"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, numba_only_t, "o-", label="NUMBA ONLY")
    plt.plot(num_paths_list, numba_parallel_t, "o-", label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, numba_only_v, "o-", label="NUMBA ONLY")
    plt.plot(num_paths_list, numba_parallel_v, "o-", label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)


########################################################################################

if 1 == 0:
    test_fin_numba_numpy_speed(False)
    test_fin_numba_numpy_speed(True)

########################################################################################

if 1 == 0:
    test_fin_numba_numba_parallel(False)
    test_fin_numba_numba_parallel(True)

########################################################################################

#  test_cases.compare_test_cases()
