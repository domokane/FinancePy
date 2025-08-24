# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys
from financepy.models.process_simulator import FinVasicekNumericalScheme
from financepy.models.process_simulator import CIRNumericalScheme
from financepy.models.process_simulator import FinHestonNumericalScheme
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import FinProcessSimulator
from FinTestCases import FinTestCases, global_test_case_mode

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_process_simulator():

    import time

    num_paths = 20000
    num_annual_steps = 100
    seed = 1919
    t = 1.0
    model_sim = FinProcessSimulator()
    print_paths = False

    test_cases.banner(
        "######################## GBM NORMAL ###############################"
    )
    sigma = 0.10
    stock_price = 100.0
    drift = 0.04
    scheme = FinGBMNumericalScheme.NORMAL_SCHEME
    model_params = (stock_price, drift, sigma, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.GBM_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.header("PROCESS", "TIME")
    test_cases.print("GBM NORMAL", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "######################## GBM ANTITHETIC ###########################"
    )
    sigma = 0.10
    stock_price = 100.0
    drift = 0.04
    scheme = FinGBMNumericalScheme.ANTITHETIC_SCHEME
    model_params = (stock_price, drift, sigma, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.GBM_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("GBM ANTITHETIC", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "###################### HESTON euler ###############################"
    )
    stock_price = 100.0
    v0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    rho = -0.9
    scheme = FinHestonNumericalScheme.EULER_SCHEME
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON EULER", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "###################### HESTON EULERLOG ############################"
    )
    stock_price = 100.0
    v0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    rho = -0.9
    scheme = FinHestonNumericalScheme.EULERLOG_SCHEME
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON EULERLOG", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "###################### HESTON QUADEXP #############################"
    )
    stock_price = 100.0
    v0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    rho = -0.9
    scheme = FinHestonNumericalScheme.QUADEXP_SCHEME
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON QUADEXP", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "######################## VASICEK NORMAL ###########################"
    )
    r0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    scheme = FinVasicekNumericalScheme.NORMAL
    model_params = (r0, kappa, theta, sigma, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.VASICEK_PROCESS,
        t,
        model_params,
        num_annual_steps,
        num_paths,
        seed,
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("VASICEK_NORMAL", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "####################### VASICEK ANTITHETIC ########################"
    )
    r0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    scheme = FinVasicekNumericalScheme.ANTITHETIC
    model_params = (r0, kappa, theta, sigma, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.VASICEK_PROCESS,
        t,
        model_params,
        num_annual_steps,
        num_paths,
        seed,
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("VASICEK_NORMAL ANTI", elapsed)
    if print_paths:
        print(paths)

    test_cases.banner(
        "############################# CIR #################################"
    )
    r0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    scheme = CIRNumericalScheme.MILSTEIN_SCHEME
    model_params = (r0, kappa, theta, sigma, scheme)
    start = time.time()
    paths = model_sim.get_process(
        ProcessTypes.CIR_PROCESS, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("CIR", elapsed)
    if print_paths:
        print(paths)


########################################################################################

test_fin_process_simulator()
test_cases.compare_test_cases()
