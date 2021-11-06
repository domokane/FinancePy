###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.stats import mean, stdev, correlation
import numpy as np
seed = 1972
np.random.seed(seed)
num_trials = 1000000
x = np.random.normal(0.0, 1.0, size=(num_trials))
y = np.random.normal(0.0, 1.0, size=(num_trials))


def test_mean():
    np_result = x.mean()
    fp_result = mean(x)
    assert round(fp_result, 10) == round(np_result, 10)


def test_stdev():
    np_result = x.std()
    fp_result = stdev(x)
    assert round(fp_result, 10) == round(np_result, 10)

# TODO: tests for stderr, var, moment


def test_correlation():
    beta = 0.4
    z = x * beta + y * np.sqrt(1.0 - beta * beta)

    np_result = np.corrcoef(x, z)[0, 1]
    fp_result = correlation(x, z)
    assert round(fp_result, 10) == round(np_result, 10)
