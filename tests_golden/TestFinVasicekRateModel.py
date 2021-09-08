###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.vasicek_mc import zero_price, zero_price_mc
import numpy as np
import time

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelRatesVasicek():

    r0 = 0.05
    a = 0.10
    b = 0.05
    sigma = 0.05
    t = 5.0

    p = zero_price(r0, a, b, sigma, t)

    num_paths = 1000
    dt = 0.02
    seed = 1968

    testCases.header("TIME", "T", "P", "P_MC", "P_MC2")

    for t in np.linspace(0, 10, 21):
        start = time.time()
        p_MC = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
        p_MC2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
        p = zero_price(r0, a, b, sigma, t)
        end = time.time()
        elapsed = end - start
        testCases.print(elapsed, t, p, p_MC, p_MC2)


###############################################################################


test_FinModelRatesVasicek()
testCases.compareTestCases()
