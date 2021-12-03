###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.utils.date import Date
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.discount_curve import DiscountCurve
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

###############################################################################


def test_FinInterpolatedForwards():

    import matplotlib.pyplot as plt

    tValues = np.array([0.0, 3.0, 5.0, 10.0])
    rValues = np.array([0.04, 0.07, 0.08, 0.09])
    df_values = np.exp(-tValues*rValues)
    tInterpValues = np.linspace(0.0, 12.0, 49)

    curve_date = Date(1, 1, 2019)

    tDates = curve_date.add_years(tValues)
    tInterpDates = curve_date.add_years(tInterpValues)

    for interp_type in InterpTypes:

        discount_curve = DiscountCurve(
            curve_date, tDates, df_values, interp_type)
        dfInterpValues = discount_curve.df(tInterpDates)
        fwdInterpValues = discount_curve.fwd(tInterpDates)
        zeroInterpValues = discount_curve.zero_rate(tInterpDates)

        if PLOT_GRAPHS:
            plt.figure(figsize=(8, 6))
            plt.plot(tValues, df_values, 'o', color='g', label="DFS:")
            plt.plot(tInterpValues, dfInterpValues, color='r',
                     label="DF:" + str(interp_type))
            plt.legend()
            plt.figure(figsize=(8, 6))
            plt.plot(tInterpValues, fwdInterpValues, color='r',
                     label="FWD:" + str(interp_type))
            plt.plot(tInterpValues, zeroInterpValues, color='b',
                     label="ZERO:" + str(interp_type))
            plt.plot(tValues, rValues, 'o', color='g',  label="ZERO RATES")
            plt.legend()

###############################################################################


test_FinInterpolatedForwards()
testCases.compareTestCases()
