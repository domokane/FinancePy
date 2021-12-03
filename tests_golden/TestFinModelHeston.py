###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time
from financepy.models.heston import Heston, HestonNumericalScheme
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


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
    dividend_yield = 0.01
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
            hestonModel = Heston(v0, kappa, theta, sigma, rho)
            for strike_price in np.linspace(95, 105, 3):
                call_option = EquityVanillaOption(
                    expiry_date, strike_price, OptionTypes.EUROPEAN_CALL)
                value_mc_Heston = hestonModel.value_mc(
                    valuation_date,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed)
                start = time.time()
                valueGatheral = hestonModel.value_gatheral(
                    valuation_date, call_option, stock_price, interest_rate, dividend_yield)
                valueLewisRouah = hestonModel.value_lewis_rouah(
                    valuation_date, call_option, stock_price, interest_rate, dividend_yield)
                valueLewis = hestonModel.value_lewis(
                    valuation_date, call_option, stock_price, interest_rate, dividend_yield)
                valueWeber = hestonModel.value_weber(
                    valuation_date, call_option, stock_price, interest_rate, dividend_yield)
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
    dividend_yield = 0.01
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
                hestonModel = Heston(v0, kappa, theta, sigma, rho)
                call_option = EquityVanillaOption(
                    expiry_date, strike_price, OptionTypes.EUROPEAN_CALL)
                valueWeber = hestonModel.value_weber(
                    valuation_date, call_option, stock_price, interest_rate, dividend_yield)

                start = time.time()

                value_mc_EULER = hestonModel.value_mc(
                    valuation_date,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.EULER)
                value_mc_EULERLOG = hestonModel.value_mc(
                    valuation_date,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.EULERLOG)
                value_mc_QUADEXP = hestonModel.value_mc(
                    valuation_date,
                    call_option,
                    stock_price,
                    interest_rate,
                    dividend_yield,
                    num_paths,
                    num_steps,
                    seed,
                    HestonNumericalScheme.QUADEXP)

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
