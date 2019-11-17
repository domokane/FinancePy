# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np
from financepy.market.volatility.FinVolatilityCurve import FinVolatilityCurve
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.finutils.FinDate import FinDate

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinVolatilityCurve():

    valueDate = FinDate(2012, 6, 20)
    expiryDate = FinDate(2012, 12, 20)
    strikes = np.linspace(70, 130, 7)
    vols = np.array([0.23, 0.24, 0.267, 0.29, 0.31, 0.33, 0.35])
    polynomial = 5
    volCurve = FinVolatilityCurve(valueDate, expiryDate,
                                  strikes, vols, polynomial)

    interpStrikes = np.linspace(50, 150, 100)
    interpVols = volCurve.volatility(interpStrikes)

    if 1 == 0:
        import matplotlib.pyplot as plt
        plt.plot(strikes, vols, 'o', interpStrikes, interpVols)
        plt.show()


test_FinVolatilityCurve()
testCases.compareTestCases()
