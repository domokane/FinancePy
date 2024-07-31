###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.math import scale
from financepy.market.curves.discount_curve_ns import DiscountCurveNS
from financepy.utils.date import Date
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

##########################################################################


def test_FinNelsonSiegelCurve():

    tau = 2.0
    times = np.linspace(0.0, 10.0, 5)
    curve_dt = Date(6, 6, 2019)
    dates = curve_dt.add_years(times)

    curve1 = DiscountCurveNS(curve_dt, 1, 0, 0, tau)
    factor1loading = curve1.zero_rate(dates)
    curve2 = DiscountCurveNS(curve_dt, 0, 1, 0, tau)
    factor2loading = curve2.zero_rate(dates)
    curve3 = DiscountCurveNS(curve_dt, 0, 0, 1, tau)
    factor3loading = curve3.zero_rate(dates)

    test_cases.header("FACTOR LOADING ON ZERO RATES")
    test_cases.print(factor1loading)
    test_cases.print(factor2loading)
    test_cases.print(factor3loading)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(factor1loading, 1), label='beta_1')
        plt.plot(times, scale(factor2loading, 1), label='beta_2')
        plt.plot(times, scale(factor3loading, 1), label='beta_3')
        plt.ylim((0, 1.05))

        plt.title('Factor Loadings in Nelson-Siegel Model')
        plt.xlabel('Time (years)')
        plt.ylabel('Loading')
        plt.legend(loc='best')

    ###########################################################################

    test_cases.header("beta_1", "beta_2", "beta_3", "ZEROS")

    beta_1 = 0.03
    beta_2 = -0.02
    beta_3 = 0.02
    curve1 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates1 = curve1.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates1)

    beta_1 = 0.04
    beta_2 = -0.02
    beta_3 = 0.02
    curve2 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates2 = curve2.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates2)

    beta_1 = 0.05
    beta_2 = -0.02
    beta_3 = 0.02
    curve3 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates3 = curve3.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates3)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.02
    curve4 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates4 = curve4.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates4)

    beta_1 = 0.07
    beta_2 = -0.02
    beta_3 = 0.02
    curve5 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates5 = curve5.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zero_rates1, 100), label='beta_1=3%')
        plt.plot(times, scale(zero_rates2, 100), label='beta_1=4%')
        plt.plot(times, scale(zero_rates3, 100), label='beta_1=5%')
        plt.plot(times, scale(zero_rates4, 100), label='beta_1=6%')
        plt.plot(times, scale(zero_rates5, 100), label='beta_1=7%')
        plt.ylim((0, 7.5))

        plt.title('Nelson-Siegel Zero Rate Curves')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    ###########################################################################

    beta_1 = 0.06
    beta_2 = -0.04
    beta_3 = 0.02
    curve1 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates1 = curve1.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates1)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.02
    curve2 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates2 = curve2.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates2)

    beta_1 = 0.06
    beta_2 = 0.00
    beta_3 = 0.02
    curve3 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates3 = curve3.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates3)

    beta_1 = 0.06
    beta_2 = 0.02
    beta_3 = 0.02
    curve4 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates4 = curve4.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates4)

    beta_1 = 0.06
    beta_2 = 0.04
    beta_3 = 0.02
    curve5 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates5 = curve5.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, zero_rates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zero_rates1, 100), label='beta_2=-4%')
        plt.plot(times, scale(zero_rates2, 100), label='beta_2=-2%')
        plt.plot(times, scale(zero_rates3, 100), label='beta_2=0%')
        plt.plot(times, scale(zero_rates4, 100), label='beta_2=2%')
        plt.plot(times, scale(zero_rates5, 100), label='beta_2=4%')
        plt.ylim((0, 10))

        plt.title('Nelson-Siegel Zero Rate Curves: Varying beta_2')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = -0.02
    curve1 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates1 = curve1.zero_rate(dates)

    test_cases.print(beta_1, beta_2, beta_3, zero_rates1)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.00
    curve2 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates2 = curve2.zero_rate(dates)

    test_cases.print(beta_1, beta_2, beta_3, zero_rates2)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.02
    curve3 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates3 = curve3.zero_rate(dates)

    test_cases.print(beta_1, beta_2, beta_3, zero_rates3)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.04
    curve4 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates4 = curve4.zero_rate(dates)

    test_cases.print(beta_1, beta_2, beta_3, zero_rates4)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = 0.06
    curve5 = DiscountCurveNS(curve_dt, beta_1, beta_2, beta_3, tau)
    zero_rates5 = curve5.zero_rate(dates)

    test_cases.print(beta_1, beta_2, beta_3, zero_rates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zero_rates1, 100), label='beta_3=-2%')
        plt.plot(times, scale(zero_rates2, 100), label='beta_3=0%')
        plt.plot(times, scale(zero_rates3, 100), label='beta_3=2%')
        plt.plot(times, scale(zero_rates4, 100), label='beta_3=4%')
        plt.plot(times, scale(zero_rates5, 100), label='beta_3=6%')
        plt.ylim((3.5, 7.5))

        plt.title('Nelson-Siegel Zero Rate Curves: Varying beta_3')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

###############################################################################


test_FinNelsonSiegelCurve()
test_cases.compareTestCases()
