########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

from financepy.models.process_simulator import FinVasicekNumericalScheme
from financepy.models.process_simulator import CIRNumericalScheme
from financepy.models.process_simulator import FinHestonNumericalScheme
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import FinProcessSimulator
from FinTestCases import FinTestCases, globalTestCaseMode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinProcessSimulator():

    import time

    num_paths = 20000
    num_annual_steps = 100
    seed = 1919
    t = 1.0
    modelSim = FinProcessSimulator()
    printPaths = False

    test_cases.banner(
        "######################## GBM NORMAL ###############################"
    )
    sigma = 0.10
    stock_price = 100.0
    drift = 0.04
    scheme = FinGBMNumericalScheme.NORMAL
    model_params = (stock_price, drift, sigma, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.GBM, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.header("PROCESS", "TIME")
    test_cases.print("GBM NORMAL", elapsed)
    if printPaths:
        print(paths)

    test_cases.banner(
        "######################## GBM ANTITHETIC ###########################"
    )
    sigma = 0.10
    stock_price = 100.0
    drift = 0.04
    scheme = FinGBMNumericalScheme.ANTITHETIC
    model_params = (stock_price, drift, sigma, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.GBM, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("GBM ANTITHETIC", elapsed)
    if printPaths:
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
    scheme = FinHestonNumericalScheme.EULER
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.HESTON, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON euler", elapsed)
    if printPaths:
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
    scheme = FinHestonNumericalScheme.EULERLOG
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.HESTON, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON EULERLOG", elapsed)
    if printPaths:
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
    scheme = FinHestonNumericalScheme.QUADEXP
    model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.HESTON, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("HESTON QUADEXP", elapsed)
    if printPaths:
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
    paths = modelSim.get_process(
        ProcessTypes.VASICEK,
        t,
        model_params,
        num_annual_steps,
        num_paths,
        seed,
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("VASICEK_NORMAL", elapsed)
    if printPaths:
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
    paths = modelSim.get_process(
        ProcessTypes.VASICEK,
        t,
        model_params,
        num_annual_steps,
        num_paths,
        seed,
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("VASICEK_NORMAL ANTI", elapsed)
    if printPaths:
        print(paths)

    test_cases.banner(
        "############################# CIR #################################"
    )
    r0 = 0.05
    kappa = 0.50
    theta = 0.05
    sigma = 0.90
    scheme = CIRNumericalScheme.MILSTEIN
    model_params = (r0, kappa, theta, sigma, scheme)
    start = time.time()
    paths = modelSim.get_process(
        ProcessTypes.CIR, t, model_params, num_annual_steps, num_paths, seed
    )
    end = time.time()
    elapsed = end - start
    test_cases.print("CIR", elapsed)
    if printPaths:
        print(paths)


########################################################################################


test_FinProcessSimulator()
test_cases.compareTestCases()
