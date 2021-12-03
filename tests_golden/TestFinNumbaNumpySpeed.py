###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np
import time
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinNumbaNumpySpeed(useSobol):

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 1999

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)

    useSobolInt = int(useSobol)

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    call_option = EquityVanillaOption(expiry_date, 100.0,
                                      OptionTypes.EUROPEAN_CALL)

    value = call_option.value(valuation_date, stock_price, discount_curve,
                              dividend_yield, model)

    num_points = 20
    v_exact = [value] * num_points

    ###########################################################################
    # DO UP TO 100K AS IT IS SLOW
    ###########################################################################

    num_paths_list = np.arange(1, num_points+1, 1) * 100000

    NONUMBA_NONUMPY_v = []
    NONUMBA_NONUMPY_t = []

    print("PURE PYTHON")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_nonumba_nonumpy(valuation_date, stock_price, discount_curve,
                                                        dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NONUMBA_NONUMPY_v.append(value_mc)
        NONUMBA_NONUMPY_t.append(duration+1e-10)

    NUMPY_ONLY_v = []
    NUMPY_ONLY_t = []

    print("NUMPY ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_only(valuation_date, stock_price, discount_curve,
                                                   dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMPY_ONLY_v.append(value_mc)
        NUMPY_ONLY_t.append(duration+1e-10)

#    speedUp = np.array(NONUMBA_NONUMPY_t)/np.array(NUMPY_ONLY_t)
#    print(NUMPY_ONLY_t)
#    print(NONUMBA_NONUMPY_t)
#    print(speedUp)

    if useSobol:
        title = "SOBOL: PURE PYTHON VS NUMPY"
    else:
        title = "PSEUDORANDOM: PURE PYTHON VS NUMPY"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, NONUMBA_NONUMPY_t, 'o-', label="PURE PYTHON")
    plt.plot(num_paths_list, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, NONUMBA_NONUMPY_v, 'o-', label="PURE PYTHON")
    plt.plot(num_paths_list, NUMPY_ONLY_v, 'o-', label="NUMPY ONLY")

    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)

    ###########################################################################
    # DO UP TO 10 MILLION NOW THAT WE HAVE NUMPY
    ###########################################################################

    num_paths_list = np.arange(1, num_points+1, 1) * 1000000

    NUMPY_ONLY_v = []
    NUMPY_ONLY_t = []

    print("NUMPY ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_only(valuation_date, stock_price, discount_curve,
                                                   dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMPY_ONLY_v.append(value_mc)
        NUMPY_ONLY_t.append(duration)

    NUMBA_NUMPY_v = []
    NUMBA_NUMPY_t = []

    print("NUMBA+NUMPY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numpy_numba(valuation_date, stock_price, discount_curve,
                                                    dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMBA_NUMPY_v.append(value_mc)
        NUMBA_NUMPY_t.append(duration)

    NUMBA_ONLY_v = []
    NUMBA_ONLY_t = []

    print("NUMBA ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_only(valuation_date, stock_price, discount_curve,
                                                   dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMBA_ONLY_v.append(value_mc)
        NUMBA_ONLY_t.append(duration)

    NUMBA_PARALLEL_v = []
    NUMBA_PARALLEL_t = []

    print("NUMBA PARALLEL")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_parallel(valuation_date, stock_price, discount_curve,
                                                       dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMBA_PARALLEL_v.append(value_mc)
        NUMBA_PARALLEL_t.append(duration)

#    speedUp = np.array(NUMBA_ONLY_t)/np.array(NUMBA_PARALLEL_t)
#    print("PARALLEL:", speedUp)

    ###########################################################################
    # COMPUTED USING NUMSTEPS FROM 1M to 10M
    ###########################################################################

    CPP_t = np.array([0.075, 0.155, 0.223, 0.313, 0.359, 0.421, 0.495, 0.556, 0.64, 0.702,
                      0.765, 0.841, 0.923, 0.982, 1.05, 1.125, 1.195, 1.261, 1.333, 1.408])

    CPP_v = np.array([9.30872, 9.29576, 9.29422, 9.29832, 9.29863, 9.30153, 9.2994, 9.3025, 9.29653, 9.29875,
                      9.29897, 9.29996, 9.29931, 9.29796, 9.29784, 9.2992, 9.3001, 9.30093, 9.29876, 9.29921])

    if useSobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.plot(num_paths_list, NUMBA_NUMPY_t, 'o-', label="NUMBA + NUMPY")

    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    ###########################################################################

    if useSobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.plot(num_paths_list, NUMBA_NUMPY_t, 'o-', label="NUMBA + NUMPY")
    plt.plot(num_paths_list, NUMBA_ONLY_t, 'o-', label="NUMBA ONLY")
    plt.plot(num_paths_list, NUMBA_PARALLEL_t, 'o-', label="NUMBA PARALLEL")

    if useSobol == False:
        plt.plot(num_paths_list, CPP_t, 'o-', label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    ###########################################################################

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, NUMBA_ONLY_v, 'o-', label="NUMBA ONLY")
    plt.plot(num_paths_list, CPP_v, 'o-', label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)

###############################################################################


def test_FinNumbaNumbaParallel(useSobol):

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 2021

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)

    useSobolInt = int(useSobol)

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    call_option = EquityVanillaOption(expiry_date, 100.0,
                                      OptionTypes.EUROPEAN_CALL)

    value = call_option.value(valuation_date, stock_price, discount_curve,
                              dividend_yield, model)

    num_points = 20
    v_exact = [value] * num_points

    num_paths_list = np.arange(1, num_points+1, 1) * 1000000

    NUMBA_ONLY_v = []
    NUMBA_ONLY_t = []

    print("NUMBA ONLY")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_only(valuation_date, stock_price, discount_curve,
                                                   dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMBA_ONLY_v.append(value_mc)
        NUMBA_ONLY_t.append(duration)

    NUMBA_PARALLEL_v = []
    NUMBA_PARALLEL_t = []

    print("NUMBA PARALLEL")
    for num_paths in num_paths_list:

        start = time.time()
        value_mc = call_option.value_mc_numba_parallel(valuation_date, stock_price, discount_curve,
                                                       dividend_yield, model, num_paths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" %
              (num_paths, value, value_mc, duration))

        NUMBA_PARALLEL_v.append(value_mc)
        NUMBA_PARALLEL_t.append(duration)

    ###########################################################################

    import matplotlib.pyplot as plt

    if useSobol:
        title = "SOBOL: NUMBA VS NUMBA + PARALLEL"
    else:
        title = "PSEUDORANDOM: NUMBA VS NUMBA + PARALLEL"

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, NUMBA_ONLY_t, 'o-', label="NUMBA ONLY")
    plt.plot(num_paths_list, NUMBA_PARALLEL_t, 'o-', label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8, 6))
    plt.plot(num_paths_list, v_exact, label="EXACT")
    plt.plot(num_paths_list, NUMBA_ONLY_v, 'o-', label="NUMBA ONLY")
    plt.plot(num_paths_list, NUMBA_PARALLEL_v, 'o-', label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)

###############################################################################


if 1 == 0:
    test_FinNumbaNumpySpeed(False)
    test_FinNumbaNumpySpeed(True)

if 1 == 0:
    test_FinNumbaNumbaParallel(False)
    test_FinNumbaNumbaParallel(True)

#  testCases.compareTestCases()
