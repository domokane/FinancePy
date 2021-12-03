###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
import matplotlib.pyplot as plt
from financepy.models.gauss_copula_onefactor import loss_dbn_recursion_gcd
from financepy.models.gauss_copula_onefactor import loss_dbn_hetero_adj_binomial
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False

##########################################################################


def test_FinLossDbnBuilder():

    num_credits = 125

    x = np.linspace(0, num_credits, num_credits + 1)
    defaultProb = 0.30
    num_steps = 25
    lossUnits = np.ones(num_credits)
    loss_ratio = np.ones(num_credits)

    testCases.header(
        "BETA",
        "BUILDER",
        "LOSS0",
        "LOSS1",
        "LOSS2",
        "LOSS3",
        "TIME")

    for beta in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:

        default_probs = np.ones(num_credits) * defaultProb
        beta_vector = np.ones(num_credits) * beta

        start = time.time()

        dbn1 = loss_dbn_recursion_gcd(num_credits,
                                      default_probs,
                                      lossUnits,
                                      beta_vector,
                                      num_steps)

        end = time.time()

        testCases.print(
            beta,
            "FULL_GCD",
            dbn1[0],
            dbn1[1],
            dbn1[2],
            dbn1[3],
            end - start)

        #######################################################################

        start = time.time()
        dbn2 = loss_dbn_hetero_adj_binomial(num_credits,
                                            default_probs,
                                            loss_ratio,
                                            beta_vector,
                                            num_steps)
        end = time.time()

        testCases.print(
            beta,
            "ADJ_BIN",
            dbn2[0],
            dbn2[1],
            dbn2[2],
            dbn2[3],
            end - start)

        #######################################################################

        if plotGraphs:
            plt.figure()
            plt.plot(x, dbn1, label='GCD FULL')
            plt.plot(x, dbn2, label='ADJ BIN')
            plt.legend()
            plt.show()

            dbn3 = dbn2 - dbn1
            plt.plot(x, dbn3, label='DIFF')
            plt.legend()
            plt.show()

    #######################################################################
    # INHOMOGENEOUS CASE
    #######################################################################

    num_credits = 100
    beta = 0.0
    defaultProb = 0.10

    default_probs = np.random.randint(3, 4, size=(num_credits)) / 10.0
    beta_vector = np.random.randint(5, 6, size=(num_credits)) / 10.0
    lossUnits = np.random.randint(1, 3, size=(num_credits)) / 1.0

    start = time.time()
    dbn1 = loss_dbn_recursion_gcd(num_credits,
                                  default_probs,
                                  lossUnits,
                                  beta_vector,
                                  num_steps)
    end = time.time()

    testCases.print(
        beta,
        "ADJ_BIN",
        dbn1[0],
        dbn1[1],
        dbn1[2],
        dbn1[3],
        end - start)

    start = time.time()
    dbn2 = loss_dbn_hetero_adj_binomial(num_credits,
                                        default_probs,
                                        loss_ratio,
                                        beta_vector,
                                        num_steps)
    end = time.time()

    testCases.print(
        beta,
        "ADJ_BIN",
        dbn2[0],
        dbn2[1],
        dbn2[2],
        dbn2[3],
        end - start)

##########################################################################


test_FinLossDbnBuilder()
testCases.compareTestCases()
