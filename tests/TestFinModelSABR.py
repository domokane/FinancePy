###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.models.FinModelSABR import blackVolFromSABR

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

alpha = 0.28
beta = 1.0
rho = -0.09
nu = 0.21

f = 0.043
k = 0.050
t = 1.0

###############################################################################


def test_SABR():

    testCases.header("ALPHA", "BETA", "RHO", "VOL")

    for alpha in [0.1, 0.2, 0.3]:
        for beta in [0.5, 1.0, 2.0]:
            for rho in [-0.8, 0.0, 0.8]:
                vol = blackVolFromSABR(alpha, beta, rho, nu, f, k, t)
                testCases.print(alpha, beta, rho, vol)

###############################################################################


test_SABR()
testCases.compareTestCases()
