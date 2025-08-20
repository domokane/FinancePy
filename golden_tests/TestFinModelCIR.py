########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys

sys.path.append("..")

import time
import numpy as np
from financepy.models.cir_montecarlo import zero_price_mc, zero_price
from financepy.models.cir_montecarlo import CIRNumericalScheme
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

########################################################################################


def test_FinModelRatesCIR():

    r0 = 0.05
    a = 0.20
    b = 0.05
    sigma = 0.20
    t = 5.0

    num_paths = 2000
    dt = 0.05
    seed = 1968

    test_cases.header(
        "MATURITY",
        "TIME",
        "FORMULA",
        "EULER",
        "LOGNORM",
        "MILSTEIN",
        "KJ",
        "EXACT",
    )

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
            CIRNumericalScheme.EULER.value,
        )
        p_MC2 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.LOGNORMAL.value,
        )
        p_MC3 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.MILSTEIN.value,
        )
        p_MC4 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.KAHLJACKEL.value,
        )
        p_MC5 = zero_price_mc(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            num_paths,
            seed,
            CIRNumericalScheme.EXACT.value,
        )
        end = time.time()
        elapsed = end - start
        test_cases.print(t, elapsed, p, p_MC1, p_MC2, p_MC3, p_MC4, p_MC5)


########################################################################################


test_FinModelRatesCIR()
test_cases.compareTestCases()
