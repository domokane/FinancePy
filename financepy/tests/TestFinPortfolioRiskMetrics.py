# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 17:04:57 2019

@author: Dominic
"""

import numpy as np

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.risk.FinPortfolioRiskMetrics import expectedLoss
from financepy.risk.FinPortfolioRiskMetrics import expectedShortfall
from financepy.risk.FinPortfolioRiskMetrics import valueAtRisk

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinPortfolioRiskMetrics():

    losses = np.array([2.0, 11.0, 15.0])
    probabilities = np.array([0.35, 0.40, 0.25])
    confidence = 0.95
    el = expectedLoss(losses, probabilities)
    var = valueAtRisk(losses, probabilities, confidence)
    es = expectedShortfall(losses, probabilities, confidence)

    print(el, var, es)


test_FinPortfolioRiskMetrics()
