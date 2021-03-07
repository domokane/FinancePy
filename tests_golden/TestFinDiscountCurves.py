###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")

from financepy.utils.date import *
from financepy.market.discount.interpolator import FinInterpTypes
from financepy.market.discount.curve import DiscountCurve
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.market.discount.curve_ns import DiscountCurveNS
from financepy.market.discount.curve_nss import DiscountCurveNSS
from financepy.market.discount.curve_pwf import DiscountCurvePWF
from financepy.market.discount.curve_pwl import DiscountCurvePWL
from financepy.market.discount.curve_zeros import DiscountCurveZeros
from financepy.market.discount.curve_poly import DiscountCurvePoly
from financepy.utils.global_vars import gDaysInYear

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

setDateFormatType(DateFormatTypes.UK_LONG)

PLOT_GRAPHS = False

###############################################################################
# TODO: Add other discount discount
###############################################################################


def test_FinDiscountCurves():

    # Create a curve from times and discount factors
    valuation_date = Date(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0]
    dates = valuation_date.addYears(years)
    years2 = []

    for dt in dates:
        y = (dt - valuation_date) / gDaysInYear
        years2.append(y)

    rates = np.array([0.05, 0.06, 0.065, 0.07, 0.075])
    discount_factors = np.exp(-np.array(rates) * np.array(years2))
    curvesList = []

    finDiscountCurve = DiscountCurve(valuation_date, dates, discount_factors,
                                     FinInterpTypes.FLAT_FWD_RATES)
    curvesList.append(finDiscountCurve)

    finDiscountCurveFlat = DiscountCurveFlat(valuation_date, 0.05)
    curvesList.append(finDiscountCurveFlat)

    finDiscountCurveNS = DiscountCurveNS(valuation_date, 0.0305, -0.01,
                                         0.08, 10.0)
    curvesList.append(finDiscountCurveNS)

    finDiscountCurveNSS = DiscountCurveNSS(valuation_date, 0.035, -0.02,
                                           0.09, 0.1, 1.0, 2.0)
    curvesList.append(finDiscountCurveNSS)

    finDiscountCurvePoly = DiscountCurvePoly(valuation_date, [0.05, 0.002,
                                                              -0.00005])
    curvesList.append(finDiscountCurvePoly)

    finDiscountCurvePWF = DiscountCurvePWF(valuation_date, dates, rates)
    curvesList.append(finDiscountCurvePWF)

    finDiscountCurvePWL = DiscountCurvePWL(valuation_date, dates, rates)
    curvesList.append(finDiscountCurvePWL)

    finDiscountCurveZeros = DiscountCurveZeros(valuation_date, dates, rates)
    curvesList.append(finDiscountCurveZeros)

    curveNames = []
    for curve in curvesList:
        curveNames.append(type(curve).__name__)

    testCases.banner("SINGLE CALLS NO VECTORS")
    testCases.header("CURVE", "DATE", "ZERO", "DF", "CCFWD", "MMFWD", "SWAP")

    years = np.linspace(1, 10, 10)
    fwdMaturityDates = valuation_date.addYears(years)

    testCases.banner("######################################################")
    testCases.banner("SINGLE CALLS")
    testCases.banner("######################################################")

    for name, curve in zip(curveNames, curvesList):
        for fwdMaturityDate in fwdMaturityDates:
            tenor = "3M"
            zeroRate = curve.zeroRate(fwdMaturityDate)
            fwd = curve.fwd(fwdMaturityDate)
            fwd_rate = curve.fwd_rate(fwdMaturityDate, tenor)
            swap_rate = curve.swap_rate(valuation_date, fwdMaturityDate)
            df = curve.df(fwdMaturityDate)

            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zeroRate),
                            "%8.7f" % (df),
                            "%7.6f" % (fwd),
                            "%7.6f" % (fwd_rate),
                            "%7.6f" % (swap_rate))

    # Examine vectorisation
    testCases.banner("######################################################")
    testCases.banner("VECTORISATIONS")
    testCases.banner("######################################################")

    for name, curve in zip(curveNames, curvesList):
        tenor = "3M"
        zeroRate = curve.zeroRate(fwdMaturityDates)
        fwd = curve.fwd(fwdMaturityDates)
        fwd_rate = curve.fwd_rate(fwdMaturityDates, tenor)
        swap_rate = curve.swap_rate(valuation_date, fwdMaturityDates)
        df = curve.df(fwdMaturityDates)

        for i in range(0, len(fwdMaturityDates)):
            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zeroRate[i]),
                            "%8.7f" % (df[i]),
                            "%7.6f" % (fwd[i]),
                            "%7.6f" % (fwd_rate[i]),
                            "%7.6f" % (swap_rate[i]))

    if PLOT_GRAPHS:

        years = np.linspace(0, 10, 121)
        years2 = years + 1.0
        fwdDates = valuation_date.addYears(years)
        fwdDates2 = valuation_date.addYears(years2)

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            zeroRates = curve.zeroRate(fwdDates)
            plt.plot(years, zeroRates, label=name)
        plt.legend()
        plt.title("Zero Rates")

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            fwd_rates = curve.fwd(fwdDates)
            plt.plot(years, fwd_rates, label=name)
        plt.legend()
        plt.title("CC Fwd Rates")

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            fwd_rates = curve.fwd_rate(fwdDates, fwdDates2)
            plt.plot(years, fwd_rates, label=name)
        plt.legend()
        plt.title("CC Fwd Rates")

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            dfs = curve.df(fwdDates)
            plt.plot(years, dfs, label=name)
        plt.legend()
        plt.title("Discount Factors")

###############################################################################


test_FinDiscountCurves()
testCases.compareTestCases()
