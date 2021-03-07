###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import time

import sys
sys.path.append("..")

from financepy.models.FinModelRatesCIR import zeroPrice_MC, zeroPrice
from financepy.models.FinModelRatesCIR import FinCIRNumericalScheme

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinModelRatesCIR():

    r0 = 0.05
    a = 0.20
    b = 0.05
    sigma = 0.20
    t = 5.0

    numPaths = 2000
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
        p = zeroPrice(r0, a, b, sigma, t)
        p_MC1 = zeroPrice_MC(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            numPaths,
            seed,
            FinCIRNumericalScheme.EULER.value)
        p_MC2 = zeroPrice_MC(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            numPaths,
            seed,
            FinCIRNumericalScheme.LOGNORMAL.value)
        p_MC3 = zeroPrice_MC(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            numPaths,
            seed,
            FinCIRNumericalScheme.MILSTEIN.value)
        p_MC4 = zeroPrice_MC(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            numPaths,
            seed,
            FinCIRNumericalScheme.KAHLJACKEL.value)
        p_MC5 = zeroPrice_MC(
            r0,
            a,
            b,
            sigma,
            t,
            dt,
            numPaths,
            seed,
            FinCIRNumericalScheme.EXACT.value)
        end = time.time()
        elapsed = end - start
        testCases.print(t, elapsed, p, p_MC1, p_MC2, p_MC3, p_MC4, p_MC5)

###############################################################################


test_FinModelRatesCIR()
testCases.compareTestCases()
