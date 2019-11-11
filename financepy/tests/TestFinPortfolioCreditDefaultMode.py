# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 21:57:13 2019

@author: Dominic
"""

from math import exp
import numpy as np

# TODO Set up test cases correctly

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode

from financepy.risk.FinPortfolioCreditDefaultMode import FinPortfolioCreditDefaultMode
from financepy.risk.FinPortfolioRiskMetrics import expectedLoss
from financepy.risk.FinPortfolioRiskMetrics import valueAtRisk
from financepy.risk.FinPortfolioRiskMetrics import expectedShortfall

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinPortfolioCreditDefaultMode():

    hazardRate = 0.05
    recoveryRate = 0.40
    beta = 0.50

    tmat = 5.0
    numCredits = 100
    weights = np.ones(numCredits)
    hazardRates = np.ones(numCredits) * hazardRate
    recoveryRates = np.ones(numCredits) * recoveryRate
    betas = np.ones(numCredits) * beta
    numPoints = 30

    portfolioModel = FinPortfolioCreditDefaultMode(weights)

    support, dbn = portfolioModel.lossDistribution(tmat,
                                                   hazardRates,
                                                   recoveryRates,
                                                   betas,
                                                   numPoints)

    import matplotlib.pyplot as plt
    plt.plot(support, dbn)
    plt.show()

    el1 = (1.0 - exp(-hazardRate*tmat)) * (1.0 - recoveryRate)
    el2 = expectedLoss(support, dbn)
    print("EXPECTED LOSS: THEORY %10.7f  ACTUAL %10.7f" % (el1, el2))

    print("CONFIDENCE  VALUE_AT_RISK  EXPECTED_SHORTFALL")
    for confidence in np.linspace(0.0, 0.99, 10):
        var = valueAtRisk(support, dbn, confidence)
        es = expectedShortfall(support, dbn, confidence)
        print("%10.4f %10.7f %10.7f" % (confidence, var, es))

    print("CONFIDENCE  VALUE_AT_RISK  EXPECTED_SHORTFALL")
    for confidence in np.linspace(0.99, 1.0, 101):
        var = valueAtRisk(support, dbn, confidence)
        es = expectedShortfall(support, dbn, confidence)
        print("%10.4f %10.7f %10.7f" % (confidence, var, es))

test_FinPortfolioCreditDefaultMode()
testCases.compareTestCases()
