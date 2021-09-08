###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.global_vars import gDaysInYear
from financepy.market.curves.discount_curve_poly import DiscountCurvePoly
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.discount_curve_pwl import DiscountCurvePWL
from financepy.market.curves.discount_curve_pwf import DiscountCurvePWF
from financepy.market.curves.discount_curve_nss import DiscountCurveNSS
from financepy.market.curves.discount_curve_ns import DiscountCurveNS
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date, set_date_format, DateFormatTypes
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

set_date_format(DateFormatTypes.UK_LONG)

PLOT_GRAPHS = False

###############################################################################
# TODO: Add other discount discount
###############################################################################


def test_FinDiscountCurves():

    # Create a curve from times and discount factors
    valuation_date = Date(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0]
    dates = valuation_date.add_years(years)
    years2 = []

    for dt in dates:
        y = (dt - valuation_date) / gDaysInYear
        years2.append(y)

    rates = np.array([0.05, 0.06, 0.065, 0.07, 0.075])
    discount_factors = np.exp(-np.array(rates) * np.array(years2))
    curvesList = []

    finDiscountCurve = DiscountCurve(valuation_date, dates, discount_factors,
                                     InterpTypes.FLAT_FWD_RATES)
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
    fwdMaturityDates = valuation_date.add_years(years)

    testCases.banner("######################################################")
    testCases.banner("SINGLE CALLS")
    testCases.banner("######################################################")

    for name, curve in zip(curveNames, curvesList):
        for fwdMaturityDate in fwdMaturityDates:
            tenor = "3M"
            zero_rate = curve.zero_rate(fwdMaturityDate)
            fwd = curve.fwd(fwdMaturityDate)
            fwd_rate = curve.fwd_rate(fwdMaturityDate, tenor)
            swap_rate = curve.swap_rate(valuation_date, fwdMaturityDate)
            df = curve.df(fwdMaturityDate)

            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zero_rate),
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
        zero_rate = curve.zero_rate(fwdMaturityDates)
        fwd = curve.fwd(fwdMaturityDates)
        fwd_rate = curve.fwd_rate(fwdMaturityDates, tenor)
        swap_rate = curve.swap_rate(valuation_date, fwdMaturityDates)
        df = curve.df(fwdMaturityDates)

        for i in range(0, len(fwdMaturityDates)):
            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zero_rate[i]),
                            "%8.7f" % (df[i]),
                            "%7.6f" % (fwd[i]),
                            "%7.6f" % (fwd_rate[i]),
                            "%7.6f" % (swap_rate[i]))

    if PLOT_GRAPHS:

        years = np.linspace(0, 10, 121)
        years2 = years + 1.0
        fwdDates = valuation_date.add_years(years)
        fwdDates2 = valuation_date.add_years(years2)

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            zero_rates = curve.zero_rate(fwdDates)
            plt.plot(years, zero_rates, label=name)
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
