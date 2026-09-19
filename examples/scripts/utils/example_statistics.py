# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import time
import numpy as np

import add_fp_to_path

from financepy.utils.stats import mean, stdev, correlation


########################################################################################


def test_fin_statistics():

    seed = 1972
    np.random.seed(seed)

    num_trials = 1000000
    x = np.random.normal(0.0, 1.0, size=num_trials)
    y = np.random.normal(0.0, 1.0, size=num_trials)

    # DO NUMPY TIMINGS

    start = time.time()

    print("l", "Mean", "SD")

    for l in range(0, 10):
        meanx1 = x.mean()
        sd1 = x.std()
        print(l, meanx1, sd1)

    end = time.time()
    elapsed = end - start

    start = time.time()

    print("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = np.corrcoef(x, z)[0, 1]
        print(beta, c)

    end = time.time()
    elapsed = end - start
    print("TIME")
    print(elapsed)

    # DO STATS TIMINGS

    print("l", "Mean", "SD")

    start = time.time()

    for l in range(0, 10):
        mean2 = mean(x)
        sd2 = stdev(x)
        print(l, mean2, sd2)

    end = time.time()
    elapsed = end - start
    print("TIME")
    print(elapsed)

    start = time.time()

    print("Corr", "Measured")

    for beta in np.linspace(0.0, 1.0, num=11):
        z = x * beta + y * np.sqrt(1.0 - beta * beta)
        c = correlation(x, z)
        print(beta, c)

    end = time.time()
    elapsed = end - start
    print("TIME")
    print(elapsed)


########################################################################################

test_fin_statistics()
