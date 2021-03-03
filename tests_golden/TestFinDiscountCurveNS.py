###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.market.curves.FinDiscountCurveNS import FinDiscountCurveNS
from financepy.finutils.FinMath import scale

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

PLOT_GRAPHS = False

##########################################################################


def test_FinNelsonSiegelCurve():

    tau = 2.0
    times = np.linspace(0.0, 10.0, 5)
    curveDate = FinDate(6, 6, 2019)
    dates = curveDate.addYears(times)

    curve1 = FinDiscountCurveNS(curveDate, 1, 0, 0, tau)
    factor1loading = curve1.zeroRate(dates)
    curve2 = FinDiscountCurveNS(curveDate, 0, 1, 0, tau)
    factor2loading = curve2.zeroRate(dates)
    curve3 = FinDiscountCurveNS(curveDate, 0, 0, 1, tau)
    factor3loading = curve3.zeroRate(dates)

    testCases.header("FACTOR LOADING ON ZERO RATES")
    testCases.print(factor1loading)
    testCases.print(factor2loading)
    testCases.print(factor3loading)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(factor1loading, 1), label='beta1')
        plt.plot(times, scale(factor2loading, 1), label='beta2')
        plt.plot(times, scale(factor3loading, 1), label='beta3')
        plt.ylim((0, 1.05))

        plt.title('Factor Loadings in Nelson-Siegel Model')
        plt.xlabel('Time (years)')
        plt.ylabel('Loading')
        plt.legend(loc='best')

    ###########################################################################

    testCases.header("BETA1", "BETA2", "BETA3", "ZEROS")

    beta1 = 0.03
    beta2 = -0.02
    beta3 = 0.02
    curve1 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates1 = curve1.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates1)

    beta1 = 0.04
    beta2 = -0.02
    beta3 = 0.02
    curve2 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates2 = curve2.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates2)

    beta1 = 0.05
    beta2 = -0.02
    beta3 = 0.02
    curve3 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates3 = curve3.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates3)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.02
    curve4 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates4 = curve4.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates4)

    beta1 = 0.07
    beta2 = -0.02
    beta3 = 0.02
    curve5 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates5 = curve5.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zeroRates1, 100), label='beta1=3%')
        plt.plot(times, scale(zeroRates2, 100), label='beta1=4%')
        plt.plot(times, scale(zeroRates3, 100), label='beta1=5%')
        plt.plot(times, scale(zeroRates4, 100), label='beta1=6%')
        plt.plot(times, scale(zeroRates5, 100), label='beta1=7%')
        plt.ylim((0, 7.5))

        plt.title('Nelson-Siegel Zero Rate Curves')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    ###########################################################################

    beta1 = 0.06
    beta2 = -0.04
    beta3 = 0.02
    curve1 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates1 = curve1.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates1)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.02
    curve2 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates2 = curve2.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates2)

    beta1 = 0.06
    beta2 = 0.00
    beta3 = 0.02
    curve3 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates3 = curve3.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates3)

    beta1 = 0.06
    beta2 = 0.02
    beta3 = 0.02
    curve4 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates4 = curve4.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates4)

    beta1 = 0.06
    beta2 = 0.04
    beta3 = 0.02
    curve5 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates5 = curve5.zeroRate(dates)
    testCases.print(beta1, beta2, beta3, zeroRates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zeroRates1, 100), label='beta2=-4%')
        plt.plot(times, scale(zeroRates2, 100), label='beta2=-2%')
        plt.plot(times, scale(zeroRates3, 100), label='beta2=0%')
        plt.plot(times, scale(zeroRates4, 100), label='beta2=2%')
        plt.plot(times, scale(zeroRates5, 100), label='beta2=4%')
        plt.ylim((0, 10))

        plt.title('Nelson-Siegel Zero Rate Curves: Varying beta2')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = -0.02
    curve1 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates1 = curve1.zeroRate(dates)

    testCases.print(beta1, beta2, beta3, zeroRates1)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.00
    curve2 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates2 = curve2.zeroRate(dates)

    testCases.print(beta1, beta2, beta3, zeroRates2)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.02
    curve3 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates3 = curve3.zeroRate(dates)

    testCases.print(beta1, beta2, beta3, zeroRates3)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.04
    curve4 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates4 = curve4.zeroRate(dates)

    testCases.print(beta1, beta2, beta3, zeroRates4)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = 0.06
    curve5 = FinDiscountCurveNS(curveDate, beta1, beta2, beta3, tau)
    zeroRates5 = curve5.zeroRate(dates)

    testCases.print(beta1, beta2, beta3, zeroRates5)

    if PLOT_GRAPHS:
        plt.figure(figsize=(6, 4))
        plt.plot(times, scale(zeroRates1, 100), label='beta3=-2%')
        plt.plot(times, scale(zeroRates2, 100), label='beta3=0%')
        plt.plot(times, scale(zeroRates3, 100), label='beta3=2%')
        plt.plot(times, scale(zeroRates4, 100), label='beta3=4%')
        plt.plot(times, scale(zeroRates5, 100), label='beta3=6%')
        plt.ylim((3.5, 7.5))

        plt.title('Nelson-Siegel Zero Rate Curves: Varying beta3')
        plt.xlabel('Time (years)')
        plt.ylabel('Zero Rate (%)')
        plt.legend(loc='lower right', frameon=False)

###############################################################################


test_FinNelsonSiegelCurve()
testCases.compareTestCases()
