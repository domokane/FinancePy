###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
import sys
sys.path.append("..")

from financepy.market.volatility.FinEquityVolCurve import FinEquityVolCurve
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinVolatilityCurve():

    valueDate = FinDate(20, 6, 2012)
    expiryDate = FinDate(20, 12, 2012)
    strikes = np.linspace(70, 130, 7)
    vols = np.array([0.23, 0.24, 0.267, 0.29, 0.31, 0.33, 0.35])
    polynomial = 5
    volCurve = FinEquityVolCurve(valueDate, expiryDate,
                                 strikes, vols, polynomial)

    interpStrikes = np.linspace(50, 150, 10)
    interpVols = volCurve.volatility(interpStrikes)

    if PLOT_GRAPHS:
        import matplotlib.pyplot as plt
        plt.plot(strikes, vols, 'o', interpStrikes, interpVols)
        plt.show()

###############################################################################


test_FinVolatilityCurve()
testCases.compareTestCases()
