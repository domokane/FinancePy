# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import matplotlib.pyplot as plt
import numpy as np

import add_fp_to_path

from financepy.utils.math import scale
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date



# TODO: Add other discount discount

PLOT_GRAPHS = False

########################################################################################


def test_fin_discount_curve():

    # Create a curve from times and discount factors
    start_dt = Date(1, 1, 2018)
    years = np.linspace(0, 10, 6)
    rate = 0.05 + 0.005 * years - 0.0003 * years * years
    dfs = np.exp(-rate * years)
    dates = start_dt.add_years(years)

    curve = DiscountCurve(start_dt, dates, dfs, InterpTypes.FLAT_FWD_RATES)

    print("T", "DF", "ZERORATE", "CC_FWD", "MM_FWD", "SURVPROB")

    plot_years = np.linspace(0, 12, 12 * 12 + 1)[1:]
    plot_dts = start_dt.add_years(plot_years)

    # Examine dependency of curve on compounding rate
    zero_rates_a = curve.zero_rate(plot_dts, FrequencyTypes.ANNUAL)
    zero_rates_s = curve.zero_rate(plot_dts, FrequencyTypes.SEMI_ANNUAL)
    zero_rates_q = curve.zero_rate(plot_dts, FrequencyTypes.QUARTERLY)
    zero_rates_m = curve.zero_rate(plot_dts, FrequencyTypes.MONTHLY)
    zero_rates_c = curve.zero_rate(plot_dts, FrequencyTypes.CONTINUOUS)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(plot_years, scale(zero_rates_a, 100), label="A")
        plt.plot(plot_years, scale(zero_rates_s, 100), label="S")
        plt.plot(plot_years, scale(zero_rates_q, 100), label="Q")
        plt.plot(plot_years, scale(zero_rates_m, 100), label="M")
        plt.plot(plot_years, scale(zero_rates_c, 100), label="C")
        plt.ylim((5, 8))

        plt.title("Discount Curves")
        plt.xlabel("Time (years)")
        plt.ylabel("Zero Rate (%)")
        plt.legend(loc="lower right", frameon=False)

    # Examine dependency of fwd curve on the interpolation scheme

    for interp in InterpTypes:

        curve = DiscountCurve(start_dt, dates, dfs, interp)
        fwd_rates = curve.fwd_rate_inst(plot_dts)
        zero_rates = curve.zero_rate(plot_dts, FrequencyTypes.ANNUAL)
        par_rates = curve.swap_rate(start_dt, plot_dts, FrequencyTypes.ANNUAL)
        par_rates2 = curve.swap_rate_old(
            start_dt, plot_dts, FrequencyTypes.ANNUAL
        )

        if PLOT_GRAPHS:
            plt.figure(figsize=(6, 4))
            plt.plot(plot_years, scale(fwd_rates, 100), label="FWD RATES")
            plt.plot(plot_years, scale(zero_rates, 100), label="ZERO RATES")
            plt.plot(plot_years, scale(par_rates, 100), label="PAR RATES")
            plt.ylim((3.0, 8.5))

            plt.title("Forward Curves using " + str(interp))
            plt.xlabel("Time (years)")
            plt.ylabel("Fwd Rate (%)")
            plt.legend(loc="lower right", frameon=False)

    ###########################################################################

    for interp in InterpTypes:

        one_bp = 0.01  # 100bp

        curve = DiscountCurve(start_dt, dates, dfs, interp)
        fwd_rates = curve.fwd_rate_inst(plot_dts)

        bumped_curve = curve.bump_parallel(one_bp)
        fwd_rates_bumped = bumped_curve.fwd_rate_inst(plot_dts)

        if PLOT_GRAPHS:
            plt.figure(figsize=(6, 4))
            plt.plot(plot_years, scale(fwd_rates, 100), label="FWD RATES")
            plt.plot(
                plot_years,
                scale(fwd_rates_bumped, 100),
                label="FWD RATES BUMP",
            )
            plt.ylim((3.0, 8.5))

            plt.title("Forward Curves using " + str(interp))
            plt.xlabel("Time (years)")
            plt.ylabel("Fwd Rate (%)")
            plt.legend(loc="lower right", frameon=False)


###############################################################################

test_fin_discount_curve()
