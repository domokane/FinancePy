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
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def testAnalyticalModels():

    # Reference see table 4.1 of Rouah book
    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 4, 2015)
    v0 = 0.05  # initial variance of volatility
    theta = 0.05  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 0.10  # volatility of variance
    rho = -0.9  # correlation
    interestRate = 0.05
    dividendYield = 0.01
    seed = 2838

    numSteps = 100
    numPaths = 20000
    stockPrice = 100.0

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
                    expiryDate, strike_price, FinOptionTypes.EUROPEAN_CALL)
                value_mc_Heston = hestonModel.value_MC(
                    valueDate,
                    callOption,
                    stockPrice,
                    interestRate,
                    dividendYield,
                    numPaths,
                    numSteps,
                    seed)
                start = time.time()
                valueGatheral = hestonModel.value_Gatheral(
                    valueDate, callOption, stockPrice, interestRate, dividendYield)
                valueLewisRouah = hestonModel.value_Lewis_Rouah(
                    valueDate, callOption, stockPrice, interestRate, dividendYield)
                valueLewis = hestonModel.value_Lewis(
                    valueDate, callOption, stockPrice, interestRate, dividendYield)
                valueWeber = hestonModel.value_Weber(
                    valueDate, callOption, stockPrice, interestRate, dividendYield)
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
    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    v0 = 0.04  # initial variance of volatility
    theta = 0.04  # long term variance
    kappa = 2.0  # speed of variance reversion
    sigma = 1.0  # volatility of variance
    rho = -0.9  # correlation
    interestRate = 0.05
    dividendYield = 0.01
    seed = 238

    stockPrice = 100.0

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
        for numSteps in [25, 50]:
            for numPaths in [10000, 20000]:
                hestonModel = FinModelHeston(v0, kappa, theta, sigma, rho)
                callOption = EquityVanillaOption(
                    expiryDate, strike_price, FinOptionTypes.EUROPEAN_CALL)
                valueWeber = hestonModel.value_Weber(
                    valueDate, callOption, stockPrice, interestRate, dividendYield)

                start = time.time()

                value_mc_EULER = hestonModel.value_MC(
                    valueDate,
                    callOption,
                    stockPrice,
                    interestRate,
                    dividendYield,
                    numPaths,
                    numSteps,
                    seed,
                    FinHestonNumericalScheme.EULER)
                value_mc_EULERLOG = hestonModel.value_MC(
                    valueDate,
                    callOption,
                    stockPrice,
                    interestRate,
                    dividendYield,
                    numPaths,
                    numSteps,
                    seed,
                    FinHestonNumericalScheme.EULERLOG)
                value_mc_QUADEXP = hestonModel.value_MC(
                    valueDate,
                    callOption,
                    stockPrice,
                    interestRate,
                    dividendYield,
                    numPaths,
                    numSteps,
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
                                numSteps,
                                numPaths,
                                valueWeber,
                                err_EULER,
                                err_EULERLOG,
                                err_QUADEXP)

##########################################################################


testAnalyticalModels()
testMonteCarlo()
testCases.compareTestCases()
