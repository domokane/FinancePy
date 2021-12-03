###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from financepy.models.cir_mc import zero_price_mc, zero_price
from financepy.models.cir_mc import CIRNumericalScheme
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelRatesCIR():

    r0 = 0.05
    a = 0.20
    b = 0.05
    sigma = 0.20
    t = 5.0

    num_paths = 2000
    dt = 0.05
    seed = 1968

    testCases.header(
        "MATURITY",
        "TIME",
        "FORMULA",
        "EULER",
        "LOGNORM",
        "MILSTEIN",
        "KJ",
        "EXACT")

    for t in np.linspace(0, 10, 21):

        start = time.time()
        p = zero_price(r0, a, b, sigma, t)
        p_MC1 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.EULER.value)
        p_MC2 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.LOGNORMAL.value)
        p_MC3 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.MILSTEIN.value)
        p_MC4 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.KAHLJACKEL.value)
        p_MC5 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.EXACT.value)
        end = time.time()
        elapsed = end - start
        testCases.print(t, elapsed, p, p_MC1, p_MC2, p_MC3, p_MC4, p_MC5)

###############################################################################


test_FinModelRatesCIR()
testCases.compareTestCases()
