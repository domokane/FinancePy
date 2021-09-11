###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.math import scale
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
# TODO: Add other discount discount
###############################################################################

PLOT_GRAPHS = False


def test_FinDiscountCurve():

    # Create a curve from times and discount factors
    start_date = Date(1, 1, 2018)
    years = np.linspace(0, 10, 6)
    rate = 0.05 + 0.005*years - 0.0003*years*years
    dfs = np.exp(-rate * years)
    dates = start_date.add_years(years)

    curve = DiscountCurve(start_date, dates, dfs, InterpTypes.FLAT_FWD_RATES)

    testCases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    plotYears = np.linspace(0, 12, 12*12+1)[1:]
    plotDates = start_date.add_years(plotYears)

    # Examine dependency of curve on compounding rate
    zero_rates_A = curve.zero_rate(plotDates, FrequencyTypes.ANNUAL)
    zero_rates_S = curve.zero_rate(plotDates, FrequencyTypes.SEMI_ANNUAL)
    zero_rates_Q = curve.zero_rate(plotDates, FrequencyTypes.QUARTERLY)
    zero_rates_M = curve.zero_rate(plotDates, FrequencyTypes.MONTHLY)
    zero_rates_C = curve.zero_rate(plotDates, FrequencyTypes.CONTINUOUS)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(plotYears, scale(zero_rates_A, 100), label='A')
        plt.plot(plotYears, scale(zero_rates_S, 100), label='S')
        plt.plot(plotYears, scale(zero_rates_Q, 100), label='Q')
        plt.plot(plotYears, scale(zero_rates_M, 100), label='M')
        plt.plot(plotYears, scale(zero_rates_C, 100), label='C')
        plt.ylim((5, 8))

        plt.title('Discount Curves')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    # Examine dependency of fwd curve on the interpolation scheme

    for interp in InterpTypes:

        curve = DiscountCurve(start_date, dates, dfs, interp)
        fwd_rates = curve.fwd(plotDates)
        zero_rates = curve.zero_rate(plotDates, FrequencyTypes.ANNUAL)
        parRates = curve.swap_rate(
            start_date, plotDates, FrequencyTypes.ANNUAL)

        if PLOT_GRAPHS:
            plt.figure(figsize=(6, 4))
            plt.plot(plotYears, scale(fwd_rates, 100), label='FWD RATES')
            plt.plot(plotYears, scale(zero_rates, 100), label='ZERO RATES')
            plt.plot(plotYears, scale(parRates, 100), label='PAR RATES')
            plt.ylim((3.0, 8.5))

            plt.title('Forward Curves using ' + str(interp))
            plt.xlabel('Time (years)')
            plt.ylabel('Fwd Rate (%)')
            plt.legend(loc='lower right', frameon=False)

###############################################################################


test_FinDiscountCurve()
testCases.compareTestCases()
