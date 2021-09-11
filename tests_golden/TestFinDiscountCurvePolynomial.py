###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_poly import DiscountCurvePoly
import numpy as np

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################
# TODO
# Inherit from FinDiscountCurve and add df method
# Put in a convention for the rate
# Use Frequency object
##############################################################################

PLOT_GRAPHS = False


def test_FinDiscountCurvePolynomial():

    times = np.linspace(0.00, 10.0, 21)
    curve_date = Date(2, 2, 2019)
    dates = curve_date.add_years(times)
    coeffs = [0.0004, -0.0001, 0.00000010]
    curve1 = DiscountCurvePoly(curve_date, coeffs)
    zeros = curve1.zero_rate(dates)
    fwds = curve1.fwd(dates)

    if PLOT_GRAPHS:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 4))
        plt.plot(times, zeros, label="Zeros")
        plt.plot(times, fwds, label="Forwards")
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate')
        plt.legend(loc='best')


test_FinDiscountCurvePolynomial()
testCases.compareTestCases()
