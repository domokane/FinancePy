###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDate import *
from financepy.market.curves.FinInterpolator import FinInterpTypes
from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.curves.FinDiscountCurveNS import FinDiscountCurveNS
from financepy.market.curves.FinDiscountCurveNSS import FinDiscountCurveNSS
from financepy.market.curves.FinDiscountCurvePWF import FinDiscountCurvePWF
from financepy.market.curves.FinDiscountCurvePWL import FinDiscountCurvePWL
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros
from financepy.market.curves.FinDiscountCurvePoly import FinDiscountCurvePoly
from financepy.finutils.FinGlobalVariables import gDaysInYear

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

setDateFormatType(FinDateFormatTypes.UK_LONG)

PLOT_GRAPHS = False

###############################################################################
# TODO: Add other discount curves
###############################################################################


def test_FinDiscountCurves():

    # Create a curve from times and discount factors
    valuationDate = FinDate(1, 1, 2018)
    years = [1.0, 2.0, 3.0, 4.0, 5.0]
    dates = valuationDate.addYears(years)
    years2 = []

    for dt in dates:
        y = (dt - valuationDate) / gDaysInYear
        years2.append(y)

    rates = np.array([0.05, 0.06, 0.065, 0.07, 0.075])
    discountFactors = np.exp(-np.array(rates) * np.array(years2))
    curvesList = []

    finDiscountCurve = FinDiscountCurve(valuationDate, dates, discountFactors,
                                        FinInterpTypes.FLAT_FWD_RATES)
    curvesList.append(finDiscountCurve)

    finDiscountCurveFlat = FinDiscountCurveFlat(valuationDate, 0.05)
    curvesList.append(finDiscountCurveFlat)

    finDiscountCurveNS = FinDiscountCurveNS(valuationDate, 0.0305, -0.01,
                                            0.08, 10.0)
    curvesList.append(finDiscountCurveNS)

    finDiscountCurveNSS = FinDiscountCurveNSS(valuationDate, 0.035, -0.02,
                                              0.09, 0.1, 1.0, 2.0)
    curvesList.append(finDiscountCurveNSS)

    finDiscountCurvePoly = FinDiscountCurvePoly(valuationDate, [0.05, 0.002,
                                                                -0.00005])
    curvesList.append(finDiscountCurvePoly)

    finDiscountCurvePWF = FinDiscountCurvePWF(valuationDate, dates, rates)
    curvesList.append(finDiscountCurvePWF)

    finDiscountCurvePWL = FinDiscountCurvePWL(valuationDate, dates, rates)
    curvesList.append(finDiscountCurvePWL)

    finDiscountCurveZeros = FinDiscountCurveZeros(valuationDate, dates, rates)
    curvesList.append(finDiscountCurveZeros)

    curveNames = []
    for curve in curvesList:
        curveNames.append(type(curve).__name__)

    testCases.banner("SINGLE CALLS NO VECTORS")
    testCases.header("CURVE", "DATE", "ZERO", "DF", "CCFWD", "MMFWD", "SWAP")

    years = np.linspace(1, 10, 10)
    fwdMaturityDates = valuationDate.addYears(years)

    testCases.banner("######################################################")
    testCases.banner("SINGLE CALLS")
    testCases.banner("######################################################")

    for name, curve in zip(curveNames, curvesList):
        for fwdMaturityDate in fwdMaturityDates:
            tenor = "3M"
            zeroRate = curve.zeroRate(fwdMaturityDate)
            fwd = curve.fwd(fwdMaturityDate)
            fwdRate = curve.fwdRate(fwdMaturityDate, tenor)
            swapRate = curve.swapRate(valuationDate, fwdMaturityDate)
            df = curve.df(fwdMaturityDate)

            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zeroRate),
                            "%8.7f" % (df),
                            "%7.6f" % (fwd),
                            "%7.6f" % (fwdRate),
                            "%7.6f" % (swapRate))

    # Examine vectorisation
    testCases.banner("######################################################")
    testCases.banner("VECTORISATIONS")
    testCases.banner("######################################################")

    for name, curve in zip(curveNames, curvesList):
        tenor = "3M"
        zeroRate = curve.zeroRate(fwdMaturityDates)
        fwd = curve.fwd(fwdMaturityDates)
        fwdRate = curve.fwdRate(fwdMaturityDates, tenor)
        swapRate = curve.swapRate(valuationDate, fwdMaturityDates)
        df = curve.df(fwdMaturityDates)

        for i in range(0, len(fwdMaturityDates)):
            testCases.print("%-20s" % name,
                            "%-12s" % fwdMaturityDate,
                            "%7.6f" % (zeroRate[i]),
                            "%8.7f" % (df[i]),
                            "%7.6f" % (fwd[i]),
                            "%7.6f" % (fwdRate[i]),
                            "%7.6f" % (swapRate[i]))

    if PLOT_GRAPHS:

        years = np.linspace(0, 10, 121)
        years2 = years + 1.0
        fwdDates = valuationDate.addYears(years)
        fwdDates2 = valuationDate.addYears(years2)

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            zeroRates = curve.zeroRate(fwdDates)
            plt.plot(years, zeroRates, label=name)
        plt.legend()
        plt.title("Zero Rates")

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            fwdRates = curve.fwd(fwdDates)
            plt.plot(years, fwdRates, label=name)
        plt.legend()
        plt.title("CC Fwd Rates")

        plt.figure()
        for name, curve in zip(curveNames, curvesList):
            fwdRates = curve.fwdRate(fwdDates, fwdDates2)
            plt.plot(years, fwdRates, label=name)
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
