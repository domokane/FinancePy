###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from financepy.utils.stats import mean, stdev, correlation
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

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

    testCases.header("l", "Mean", "SD")

    for l in range(0, 10):
        meanx1 = x.mean()
        sd1 = x.std()
        testCases.print(l, meanx1, sd1)

    end = time.time()
    elapsed = end - start

    start = time.time()

    testCases.header("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = np.corrcoef(x, z)[0, 1]
        testCases.print(beta, c)

    end = time.time()
    elapsed = end - start
    testCases.header("TIME")
    testCases.print(elapsed)

    ##########################################################################
    # DO STATS TIMINGS

    testCases.header("l", "Mean", "SD")

    start = time.time()

    for l in range(0, 10):
        mean2 = mean(x)
        sd2 = stdev(x)
        testCases.print(l, mean2, sd2)

    end = time.time()
    elapsed = end - start
    testCases.header("TIME")
    testCases.print(elapsed)

    start = time.time()

    testCases.header("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = correlation(x, z)
        testCases.print(beta, c)

    end = time.time()
    elapsed = end - start
    testCases.header("TIME")
    testCases.print(elapsed)

###############################################################################


test_FinStatistics()
testCases.compareTestCases()
