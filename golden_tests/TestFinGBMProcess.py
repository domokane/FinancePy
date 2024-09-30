###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from os.path import dirname, join

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_types import SwapTypes
from financepy.models.gbm_process_simulator import get_assets_paths_times
from financepy.utils.math import corr_matrix_generator

test_cases = FinTestCases(__file__, globalTestCaseMode)


def testFinGBMProcess():

    num_assets = 3
    num_paths = 6
    num_time_steps = 1
    t = 1.0
    mus = 0.03 * np.ones(num_assets)
    stock_prices = 100.0 * np.ones(num_assets)
    volatilities = 0.2 * np.ones(num_assets)
    rho = 0.8
    corr_matrix = corr_matrix_generator(rho, num_assets)
    seed = 1912

    times, paths = get_assets_paths_times(
        num_assets,
        num_paths,
        num_time_steps,
        t,
        mus,
        stock_prices,
        volatilities,
        corr_matrix,
        seed,
    )


###############################################################################


testFinGBMProcess()
test_cases.compareTestCases()
