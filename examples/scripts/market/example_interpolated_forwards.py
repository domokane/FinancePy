# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import numpy as np

import add_fp_to_path

from financepy.utils.date import Date
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.discount_curve import DiscountCurve



PLOT_GRAPHS = False

########################################################################################


def test_fin_interpolated_forwards():

    import matplotlib.pyplot as plt

    t_values = np.array([0.0, 3.0, 5.0, 10.0])
    r_values = np.array([0.04, 0.07, 0.08, 0.09])
    df_values = np.exp(-t_values * r_values)
    t_interp_values = np.linspace(0.0, 12.0, 49)

    curve_dt = Date(1, 1, 2019)

    t_dates = curve_dt.add_years(t_values)
    t_interp_dates = curve_dt.add_years(t_interp_values)

    for interp_type in InterpTypes:

        discount_curve = DiscountCurve(
            curve_dt, t_dates, df_values, interp_type
        )
        df_interp_values = discount_curve.df(t_interp_dates)
        fwd_interp_values = discount_curve.fwd_rate_inst(t_interp_dates)
        zero_interp_values = discount_curve.zero_rate(t_interp_dates)

        if PLOT_GRAPHS:
            plt.figure(figsize=(8, 6))
            plt.plot(t_values, df_values, "o", color="g", label="DFS:")
            plt.plot(
                t_interp_values,
                df_interp_values,
                color="r",
                label="DF:" + str(interp_type),
            )
            plt.legend()
            plt.figure(figsize=(8, 6))
            plt.plot(
                t_interp_values,
                fwd_interp_values,
                color="r",
                label="FWD:" + str(interp_type),
            )
            plt.plot(
                t_interp_values,
                zero_interp_values,
                color="b",
                label="ZERO:" + str(interp_type),
            )
            plt.plot(t_values, r_values, "o", color="g", label="ZERO RATES")
            plt.legend()


########################################################################################

test_fin_interpolated_forwards()
