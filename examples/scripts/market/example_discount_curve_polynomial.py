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
import matplotlib.pyplot as plt

import add_fp_to_path

from financepy.utils.date import Date
from financepy.market.curves.poly_discount_curve import PolyDiscountCurve



# TODO
# Inherit from DiscountCurve and add df method
# Put in a convention for the rate
# Use Frequency object

PLOT_GRAPHS = False

########################################################################################


def test_fin_discount_curve_polynomial():

    times = np.linspace(0.00, 10.0, 21)
    curve_dt = Date(2, 2, 2019)
    dates = curve_dt.add_years(times)
    coeffs = [0.0004, -0.0001, 0.00000010]
    curve1 = PolyDiscountCurve(curve_dt, coeffs)
    zeros = curve1.zero_rate(dates)
    fwds = curve1.fwd_rate_inst(dates)

    if PLOT_GRAPHS:

        plt.figure(figsize=(6, 4))
        plt.plot(times, zeros, label="Zeros")
        plt.plot(times, fwds, label="Forwards")
        plt.xlabel("Time (years)")
        plt.ylabel("Zero Rate")
        plt.legend(loc="best")


########################################################################################

test_fin_discount_curve_polynomial()
