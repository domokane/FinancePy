########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys
sys.path.append("..")

import matplotlib.pyplot as plt
import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.math import scale
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date

test_cases = FinTestCases(__file__, globalTestCaseMode)

########################################################################################
# TODO: Add other discount discount
########################################################################################

plot_graphs = False


def test_FinDiscountCurve():

    # Create a curve from times and discount factors
    start_dt = Date(1, 1, 2018)
    years = np.linspace(0, 10, 6)
    rate = 0.05 + 0.005 * years - 0.0003 * years * years
    dfs = np.exp(-rate * years)
    dates = start_dt.add_years(years)

    curve = DiscountCurve(start_dt, dates, dfs, InterpTypes.FLAT_FWD_RATES)

    test_cases.header("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    plot_years = np.linspace(0, 12, 12 * 12 + 1)[1:]
    plot_dts = start_dt.add_years(plot_years)

    # Examine dependency of curve on compounding rate
    zero_rates_A = curve.zero_rate(plot_dts, FrequencyTypes.ANNUAL)
    zero_rates_S = curve.zero_rate(plot_dts, FrequencyTypes.SEMI_ANNUAL)
    zero_rates_Q = curve.zero_rate(plot_dts, FrequencyTypes.QUARTERLY)
    zero_rates_M = curve.zero_rate(plot_dts, FrequencyTypes.MONTHLY)
    zero_rates_C = curve.zero_rate(plot_dts, FrequencyTypes.CONTINUOUS)

    if plot_graphs:
        plt.figure(figsize=(6, 4))
        plt.plot(plot_years, scale(zero_rates_A, 100), label="A")
        plt.plot(plot_years, scale(zero_rates_S, 100), label="S")
        plt.plot(plot_years, scale(zero_rates_Q, 100), label="Q")
        plt.plot(plot_years, scale(zero_rates_M, 100), label="M")
        plt.plot(plot_years, scale(zero_rates_C, 100), label="C")
        plt.ylim((5, 8))

        plt.title("Discount Curves")
        plt.xlabel("Time (years)")
        plt.ylabel("Zero Rate (%)")
        plt.legend(loc="lower right", frameon=False)

    # Examine dependency of fwd curve on the interpolation scheme

    for interp in InterpTypes:

        curve = DiscountCurve(start_dt, dates, dfs, interp)
        fwd_rates = curve.fwd(plot_dts)
        zero_rates = curve.zero_rate(plot_dts, FrequencyTypes.ANNUAL)
        par_rates = curve.swap_rate(start_dt, plot_dts, FrequencyTypes.ANNUAL)

        if plot_graphs:
            plt.figure(figsize=(6, 4))
            plt.plot(plot_years, scale(fwd_rates, 100), label="FWD RATES")
            plt.plot(plot_years, scale(zero_rates, 100), label="ZERO RATES")
            plt.plot(plot_years, scale(par_rates, 100), label="PAR RATES")
            plt.ylim((3.0, 8.5))

            plt.title("Forward Curves using " + str(interp))
            plt.xlabel("Time (years)")
            plt.ylabel("Fwd Rate (%)")
            plt.legend(loc="lower right", frameon=False)


########################################################################################


test_FinDiscountCurve()
test_cases.compareTestCases()
