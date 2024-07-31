###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.math import scale
from financepy.market.curves.discount_curve_nss import DiscountCurveNSS
from financepy.utils.date import Date
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

##########################################################################


def test_FinNelsonSiegelSvenssonCurve():

    tau_1 = 2.0
    tau_2 = 0.5
    times = np.linspace(0.0, 10.0, 5)
    start_dt = Date(1, 1, 2020)
    dates = start_dt.add_years(times)

    curve1 = DiscountCurveNSS(start_dt, 1., 0., 0., 0., tau_1, tau_2)
    factor1loading = curve1.zero_rate(dates)
    curve2 = DiscountCurveNSS(start_dt, 0., 1., 0., 0., tau_1, tau_2)
    factor2loading = curve2.zero_rate(dates)
    curve3 = DiscountCurveNSS(start_dt, 0., 0., 1., 0., tau_1, tau_2)
    factor3loading = curve3.zero_rate(dates)
    curve4 = DiscountCurveNSS(start_dt, 0., 0., 0., 1., tau_1, tau_2)
    factor4loading = curve4.zero_rate(dates)

    test_cases.header("FACTOR LOADING ON ZERO RATES")
    test_cases.print(factor1loading)
    test_cases.print(factor2loading)
    test_cases.print(factor3loading)
    test_cases.print(factor4loading)

#    plt.figure(figsize = (6,4))
#    plt.plot(times,scaleVector(factor1loading,1),label='beta_1');
#    plt.plot(times,scaleVector(factor2loading,1),label='beta_2');
#    plt.plot(times,scaleVector(factor3loading,1),label='beta_3');
#    plt.ylim((0,1.05))
#
#    plt.title('Factor Loadings in Nelson-Siegel Model');
#    plt.xlabel('Time (years)');
#    plt.ylabel('Loading');
#    plt.legend(loc='best')

##########################################################################

    test_cases.header("beta_1", "beta_2", "beta_3", "beta_4", "ZEROS")

    beta_1 = 0.03
    beta_2 = -0.02
    beta_3 = -0.02
    beta_4 = 0.08
    curve1 = DiscountCurveNSS(start_dt,
                              beta_1, beta_2, beta_3, beta_4, tau_1, tau_2)
    zero_rates1 = curve1.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, beta_4, zero_rates1)

    beta_1 = 0.04
    beta_2 = -0.02
    beta_3 = -0.02
    beta_4 = 0.08
    curve2 = DiscountCurveNSS(start_dt,
                              beta_1, beta_2, beta_3, beta_4, tau_1, tau_2)
    zero_rates2 = curve2.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, beta_4, zero_rates2)

    beta_1 = 0.05
    beta_2 = -0.02
    beta_3 = -0.02
    beta_4 = 0.08
    curve3 = DiscountCurveNSS(start_dt,
                              beta_1, beta_2, beta_3, beta_4, tau_1, tau_2)
    zero_rates3 = curve3.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, beta_4, zero_rates3)

    beta_1 = 0.06
    beta_2 = -0.02
    beta_3 = -0.02
    beta_4 = 0.08
    curve4 = DiscountCurveNSS(start_dt,
                              beta_1, beta_2, beta_3, beta_4, tau_1, tau_2)
    zero_rates4 = curve4.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, beta_4, zero_rates4)

    beta_1 = 0.07
    beta_2 = -0.02
    beta_3 = -0.02
    beta_4 = 0.08
    curve5 = DiscountCurveNSS(start_dt,
                              beta_1, beta_2, beta_3, beta_4, tau_1, tau_2)
    zero_rates5 = curve5.zero_rate(dates)
    test_cases.print(beta_1, beta_2, beta_3, beta_4, zero_rates5)

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

##########################################################################


test_FinNelsonSiegelSvenssonCurve()
test_cases.compareTestCases()
