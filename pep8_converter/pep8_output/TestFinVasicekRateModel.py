# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")
import time

import numpy as np
from financepy.models.vasicek_mc import zero_price, zero_price_mc
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test__fin_model_rates_vasicek():

    r0 = 0.05
    a = 0.10
    b = 0.05
    sigma = 0.05
    t = 5.0

    p = zero_price(r0, a, b, sigma, t)

    num_paths = 1000
    dt = 0.02
    seed = 1968

    test_cases.header("TIME", "T", "P", "P_MC", "P_MC2")

    for t in np.linspace(0, 10, 21):
        start = time.time()
        p_mc = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
        p_mc2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
        p = zero_price(r0, a, b, sigma, t)
        end = time.time()
        elapsed = end - start
        test_cases.print(elapsed, t, p, p_mc, p_mc2)




################################################################################

test__fin_model_rates_vasicek()
test_cases.compare_test_cases()
