###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinVolatilityCurve():

    value_dt = Date(20, 6, 2012)
    expiry_dt = Date(20, 12, 2012)
    strikes = np.linspace(70, 130, 7)
    vols = np.array([0.23, 0.24, 0.267, 0.29, 0.31, 0.33, 0.35])
    polynomial = 5
    vol_curve = EquityVolCurve(value_dt, expiry_dt, strikes, vols, polynomial)

    interp_strikes = np.linspace(50, 150, 10)
    interp_vols = vol_curve.volatility(interp_strikes)

    if PLOT_GRAPHS:
        import matplotlib.pyplot as plt

        plt.plot(strikes, vols, "o", interp_strikes, interp_vols)
        plt.show()


###############################################################################


test_FinVolatilityCurve()
test_cases.compareTestCases()
