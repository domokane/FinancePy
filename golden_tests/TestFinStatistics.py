###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from financepy.utils.stats import mean, stdev, correlation
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinStatistics():
    seed = 1972
    np.random.seed(seed)

    num_trials = 1000000
    x = np.random.normal(0.0, 1.0, size=(num_trials))
    y = np.random.normal(0.0, 1.0, size=(num_trials))

    ##########################################################################
    # DO NUMPY TIMINGS
    ##########################################################################

    start = time.time()

    test_cases.header("l", "Mean", "SD")

    for l in range(0, 10):
        meanx1 = x.mean()
        sd1 = x.std()
        test_cases.print(l, meanx1, sd1)

    end = time.time()
    elapsed = end - start

    start = time.time()

    test_cases.header("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = np.corrcoef(x, z)[0, 1]
        test_cases.print(beta, c)

    end = time.time()
    elapsed = end - start
    test_cases.header("TIME")
    test_cases.print(elapsed)

    ##########################################################################
    # DO STATS TIMINGS

    test_cases.header("l", "Mean", "SD")

    start = time.time()

    for l in range(0, 10):
        mean2 = mean(x)
        sd2 = stdev(x)
        test_cases.print(l, mean2, sd2)

    end = time.time()
    elapsed = end - start
    test_cases.header("TIME")
    test_cases.print(elapsed)

    start = time.time()

    test_cases.header("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = correlation(x, z)
        test_cases.print(beta, c)

    end = time.time()
    elapsed = end - start
    test_cases.header("TIME")
    test_cases.print(elapsed)

###############################################################################


test_FinStatistics()
test_cases.compareTestCases()
