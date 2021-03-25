###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np

import sys
sys.path.append("..")

from financepy.models.heston import FinModelHeston, FinHestonNumericalScheme
from financepy.utils.global_types import FinOptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def testAnalyticalModels():

    # Reference see table 4.1 of Rouah book
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 4, 2015)
    v0 = 0.05  # initial variance of volatility
    theta = 0.05  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 0.10  # volatility of variance
    rho = -0.9  # correlation
    interest_rate = 0.05
    dividendYield = 0.01
    seed = 2838

    num_steps = 100
    num_paths = 20000
    stock_price = 100.0

    testCases.header(
        "TIME",
        "RHO",
        "SIGMA",
        "K",
        "MC",
        "GATH",
        "LEWROU",
        "LEWIS",
        "WEBER",
        "MCERR")

    for sigma in [0.5, 0.75, 1.0]:
        for rho in [-0.9, -0.5, 0.0]:
            hestonModel = FinModelHeston(v0, kappa, theta, sigma, rho)
            for strike_price in np.linspace(95, 105, 3):
                callOption = EquityVanillaOption(
                    expiry_date, strike_price, FinOptionTypes.EUROPEAN_CALL)
                value_mc_Heston = hestonModel.value_MC(
                    valuation_date,
                    callOption,
                    stock_price,
                    interest_rate,
                    dividendYield,
                    num_paths,
                    num_steps,
                    seed)
                start = time.time()
                valueGatheral = hestonModel.value_Gatheral(
                    valuation_date, callOption, stock_price, interest_rate, dividendYield)
                valueLewisRouah = hestonModel.value_Lewis_Rouah(
                    valuation_date, callOption, stock_price, interest_rate, dividendYield)
                valueLewis = hestonModel.value_Lewis(
                    valuation_date, callOption, stock_price, interest_rate, dividendYield)
                valueWeber = hestonModel.value_Weber(
                    valuation_date, callOption, stock_price, interest_rate, dividendYield)
                err = (value_mc_Heston - valueWeber)
                end = time.time()
                elapsed = end - start
                testCases.print("%6.3f" % elapsed,
                                "% 7.5f" % rho,
                                "%7.5f" % sigma,
                                "%7.2f" % strike_price,
                                "%12.9f" % value_mc_Heston,
                                "%12.9f" % valueGatheral,  # problem
                                "%12.9f" % valueLewisRouah,
                                "%12.9f" % valueLewis,
                                "%12.9f" % valueWeber,
                                "%12.9f" % err)

##########################################################################


def testMonteCarlo():

    import time

    # Reference see table 4.1 of Rouah book
    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    v0 = 0.04  # initial variance of volatility
    theta = 0.04  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 1.0  # volatility of variance
    rho = -0.9  # correlation
    interest_rate = 0.05
    dividendYield = 0.01
    seed = 238

    stock_price = 100.0

    testCases.header(
        "TIME",
        "RHO",
        "SIGMA",
        "K",
        "NSTEPS",
        "NPATHS",
        "FORMULA",
        "EULER_ERR",
        "EULLOG_ERR",
        "QE_ERR")

    for strike_price in np.linspace(95, 105, 3):
        for num_steps in [25, 50]:
            for num_paths in [10000, 20000]:
                hestonModel = FinModelHeston(v0, kappa, theta, sigma, rho)
                callOption = EquityVanillaOption(
                    expiry_date, strike_price, FinOptionTypes.EUROPEAN_CALL)
                valueWeber = hestonModel.value_Weber(
                    valuation_date, callOption, stock_price, interest_rate, dividendYield)

                start = time.time()

                value_mc_EULER = hestonModel.value_MC(
                    valuation_date,
                    callOption,
                    stock_price,
                    interest_rate,
                    dividendYield,
                    num_paths,
                    num_steps,
                    seed,
                    FinHestonNumericalScheme.EULER)
                value_mc_EULERLOG = hestonModel.value_MC(
                    valuation_date,
                    callOption,
                    stock_price,
                    interest_rate,
                    dividendYield,
                    num_paths,
                    num_steps,
                    seed,
                    FinHestonNumericalScheme.EULERLOG)
                value_mc_QUADEXP = hestonModel.value_MC(
                    valuation_date,
                    callOption,
                    stock_price,
                    interest_rate,
                    dividendYield,
                    num_paths,
                    num_steps,
                    seed,
                    FinHestonNumericalScheme.QUADEXP)

                err_EULER = (value_mc_EULER - valueWeber)
                err_EULERLOG = (value_mc_EULERLOG - valueWeber)
                err_QUADEXP = (value_mc_QUADEXP - valueWeber)

                end = time.time()
                elapsed = end - start

                testCases.print(elapsed, rho,
                                sigma,
                                strike_price,
                                num_steps,
                                num_paths,
                                valueWeber,
                                err_EULER,
                                err_EULERLOG,
                                err_QUADEXP)

##########################################################################


testAnalyticalModels()
testMonteCarlo()
testCases.compareTestCases()
