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


testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

##########################################################################


def test_FinNelsonSiegelSvenssonCurve():

    tau1 = 2.0
    tau2 = 0.5
    times = np.linspace(0.0, 10.0, 5)
    start_date = Date(1, 1, 2020)
    dates = start_date.add_years(times)

    curve1 = DiscountCurveNSS(start_date, 1., 0., 0., 0., tau1, tau2)
    factor1loading = curve1.zero_rate(dates)
    curve2 = DiscountCurveNSS(start_date, 0., 1., 0., 0., tau1, tau2)
    factor2loading = curve2.zero_rate(dates)
    curve3 = DiscountCurveNSS(start_date, 0., 0., 1., 0., tau1, tau2)
    factor3loading = curve3.zero_rate(dates)
    curve4 = DiscountCurveNSS(start_date, 0., 0., 0., 1., tau1, tau2)
    factor4loading = curve4.zero_rate(dates)

    testCases.header("FACTOR LOADING ON ZERO RATES")
    testCases.print(factor1loading)
    testCases.print(factor2loading)
    testCases.print(factor3loading)
    testCases.print(factor4loading)

#    plt.figure(figsize = (6,4))
#    plt.plot(times,scaleVector(factor1loading,1),label='beta1');
#    plt.plot(times,scaleVector(factor2loading,1),label='beta2');
#    plt.plot(times,scaleVector(factor3loading,1),label='beta3');
#    plt.ylim((0,1.05))
#
#    plt.title('Factor Loadings in Nelson-Siegel Model');
#    plt.xlabel('Time (years)');
#    plt.ylabel('Loading');
#    plt.legend(loc='best')

##########################################################################

    testCases.header("BETA1", "BETA2", "BETA3", "BETA4", "ZEROS")

    beta1 = 0.03
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve1 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates1 = curve1.zero_rate(dates)
    testCases.print(beta1, beta2, beta3, beta4, zero_rates1)

    beta1 = 0.04
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve2 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates2 = curve2.zero_rate(dates)
    testCases.print(beta1, beta2, beta3, beta4, zero_rates2)

    beta1 = 0.05
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve3 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates3 = curve3.zero_rate(dates)
    testCases.print(beta1, beta2, beta3, beta4, zero_rates3)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve4 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates4 = curve4.zero_rate(dates)
    testCases.print(beta1, beta2, beta3, beta4, zero_rates4)

    beta1 = 0.07
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve5 = DiscountCurveNSS(start_date,
                              beta1, beta2, beta3, beta4, tau1, tau2)
    zero_rates5 = curve5.zero_rate(dates)
    testCases.print(beta1, beta2, beta3, beta4, zero_rates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zero_rates1, 100), label='beta1=3%')
        plt.plot(times, scale(zero_rates2, 100), label='beta1=4%')
        plt.plot(times, scale(zero_rates3, 100), label='beta1=5%')
        plt.plot(times, scale(zero_rates4, 100), label='beta1=6%')
        plt.plot(times, scale(zero_rates5, 100), label='beta1=7%')
        plt.ylim((0, 7.5))

        plt.title('Nelson-Siegel Zero Rate Curves')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

##########################################################################


test_FinNelsonSiegelSvenssonCurve()
testCases.compareTestCases()
