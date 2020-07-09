# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 17:04:57 2019

@author: Dominic
"""

import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.risk.FinPortfolioRiskMetrics import expectedLoss
from financepy.risk.FinPortfolioRiskMetrics import expectedShortfall
from financepy.risk.FinPortfolioRiskMetrics import valueAtRisk

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinPortfolioRiskMetrics():

    losses = np.array([2.0, 11.0, 15.0, 20.0])
    probabilities = np.array([0.32, 0.61, 0.04, 0.03])

    testCases.header("CONFIDENCE", "EL", "VAR", "ES")

    for confidence in np.linspace(0.01, 1.0, 100):
        el = expectedLoss(losses, probabilities)
        var = valueAtRisk(losses, probabilities, confidence)
        es = expectedShortfall(losses, probabilities, confidence)
        testCases.print(confidence, el, var, es)


test_FinPortfolioRiskMetrics()
testCases.compareTestCases()
