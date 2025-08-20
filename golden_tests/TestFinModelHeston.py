########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys

sys.path.append("..")

import numpy as np
import time
from financepy.models.heston import Heston, HestonNumericalScheme
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def testAnalyticalModels():

    # Reference see table 4.1 of Rouah book
    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 4, 2015)
    v0 = 0.05  # initial variance of volatility
    theta = 0.05  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 0.10  # volatility of variance
    rho = -0.9  # correlation
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 2838

    num_steps = 100
    num_paths = 20000
    stock_price = 100.0

    test_cases.header(
        "TIME",
        "RHO",
        "SIGMA",
        "K",
        "MC",
        "GATH",
        "LEWROU",
        "LEWIS",
        "WEBER",
        "MCERR",
    )

    for sigma in [0.5, 0.75, 1.0]:
        for rho in [-0.9, -0.5, 0.0]:
            heston_model = Heston(v0, kappa, theta, sigma, rho)
            for strike_price in np.linspace(95, 105, 3):
                call_option = EquityVanillaOption(
                    expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
                )
                value_mc_heston = heston_model.value_mc(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                )
                start = time.time()
                value_gatheral = heston_model.value_gatheral(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                )
                value_lewis_rouah = heston_model.value_lewis_rouah(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                )
                value_lewis = heston_model.value_lewis(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                )
                value_weber = heston_model.value_weber(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                )
                err = value_mc_heston - value_weber
                end = time.time()
                elapsed = end - start
                test_cases.print(
                    "%6.3f" % elapsed,
                    "% 7.5f" % rho,
                    "%7.5f" % sigma,
                    "%7.2f" % strike_price,
                    "%12.9f" % value_mc_heston,
                    "%12.9f" % value_gatheral,  # problem
                    "%12.9f" % value_lewis_rouah,
                    "%12.9f" % value_lewis,
                    "%12.9f" % value_weber,
                    "%12.9f" % err,
                )


##########################################################################


def testMonteCarlo():

    import time

    # Reference see table 4.1 of Rouah book
    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    v0 = 0.04  # initial variance of volatility
    theta = 0.04  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 1.0  # volatility of variance
    rho = -0.9  # correlation
    interest_rate = 0.05
    dividend_yield = 0.01
    seed = 238

    stock_price = 100.0

    test_cases.header(
        "TIME",
        "RHO",
        "SIGMA",
        "K",
        "NSTEPS",
        "NPATHS",
        "FORMULA",
        "EULER_ERR",
        "EULLOG_ERR",
        "QE_ERR",
    )

    for strike_price in np.linspace(95, 105, 3):
        for num_steps in [25, 50]:
            for num_paths in [10000, 20000]:
                heston_model = Heston(v0, kappa, theta, sigma, rho)
                call_option = EquityVanillaOption(
                    expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL
                )
                value_weber = heston_model.value_weber(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                )

                start = time.time()

                value_mc_euler = heston_model.value_mc(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.EULER,
                )
                value_mc_euler_log = heston_model.value_mc(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.EULERLOG,
                )
                value_mc_quadexp = heston_model.value_mc(
                    value_dt,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.QUADEXP,
                )

                err_euler = value_mc_euler - value_weber
                err_euler_log = value_mc_euler_log - value_weber
                err_quadexp = value_mc_quadexp - value_weber

                end = time.time()
                elapsed = end - start

                test_cases.print(
                    elapsed,
                    rho,
                    sigma,
                    strike_price,
                    num_steps,
                    num_paths,
                    value_weber,
                    err_euler,
                    err_euler_log,
                    err_quadexp,
                )


##########################################################################


testAnalyticalModels()
testMonteCarlo()
test_cases.compareTestCases()
