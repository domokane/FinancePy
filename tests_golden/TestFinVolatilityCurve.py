###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinVolatilityCurve():

    valuation_date = Date(20, 6, 2012)
    expiry_date = Date(20, 12, 2012)
    strikes = np.linspace(70, 130, 7)
    vols = np.array([0.23, 0.24, 0.267, 0.29, 0.31, 0.33, 0.35])
    polynomial = 5
    volCurve = EquityVolCurve(valuation_date, expiry_date,
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
