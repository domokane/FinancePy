# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.market.curves.FinInterpolate import FinInterpMethods

from financepy.market.curves.FinDiscountCurve import FinDiscountCurve
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.curves.FinDiscountCurveNS import FinDiscountCurveNS
from financepy.market.curves.FinDiscountCurveNSS import FinDiscountCurveNSS
from financepy.market.curves.FinDiscountCurvePWF import FinDiscountCurvePWF
from financepy.market.curves.FinDiscountCurvePWL import FinDiscountCurvePWL
from financepy.market.curves.FinDiscountCurveZeros import FinDiscountCurveZeros
from financepy.market.curves.FinDiscountCurvePoly import FinDiscountCurvePoly
from financepy.finutils.FinGlobalVariables import gDaysInYear

import matplotlib.pyplot as plt

import numpy as np
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

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

    finDiscountCurve = FinDiscountCurve(valuationDate, dates, discountFactors, FinInterpMethods.FLAT_FORWARDS)
    curvesList.append(finDiscountCurve)

    finDiscountCurveFlat = FinDiscountCurveFlat(valuationDate, 0.05)
    curvesList.append(finDiscountCurveFlat)

    finDiscountCurveNS = FinDiscountCurveNS(valuationDate, 0.0305, -0.01, 0.08, 10.0)
    curvesList.append(finDiscountCurveNS)

    finDiscountCurveNSS = FinDiscountCurveNSS(valuationDate, 0.035, -0.02, 0.09, 0.1, 1.0, 2.0)
    curvesList.append(finDiscountCurveNSS)

    finDiscountCurvePoly = FinDiscountCurvePoly(valuationDate, [0.05, 0.002, -0.00005])
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

    fwdDate1 = valuationDate.addYears(1)
    fwdDate2 = fwdDate1.addTenor("3M")
 
    # Do single calls (no vectorisations)

    print("")
    print("                     CURVE           1YR        2YR        3YR        4YR        5YR")

    for name, curve in zip(curveNames,curvesList):
        print(" %30s"% name, end="")
        for date in dates:
            zeroRate = curve.zeroRate(date)
            print(" %10.5f"% (zeroRate), end="")
        print("")
    
    print("")
    print("                     CURVE           ZERO      CCFWD      MMFWD       PAR         DF")

    for name, curve in zip(curveNames,curvesList):

         zeroRate = curve.zeroRate(fwdDate1)
         print(" %30s"% name, end="")

         print(" %10.5f"% (zeroRate), end="")
    
         fwd = curve.fwd(fwdDate1)
         print(" %10.5f"% (fwd), end="")

         fwdRate = curve.fwdRate(fwdDate1, fwdDate2)
         print(" %10.5f"% (fwdRate), end="")

         parRate = curve.parRate(fwdDate1)
         print(" %10.5f"% (parRate), end="")

         df = curve.df(fwdDate1)
         print(" %10.5f"% (df))

    for curve in curvesList:
        print(curve)
        
    years = np.linspace(0, 10, 121)
    years2 = years + 3
    fwdDates = valuationDate.addYears(years)
    fwdDates2 = valuationDate.addYears(years2)

    plt.figure()
    for name, curve in zip(curveNames,curvesList):
        zeroRates = curve.zeroRate(fwdDates)
        plt.plot(years, zeroRates, label=name)
    plt.legend()
    plt.title("Zero Rates")

    plt.figure()
    for name, curve in zip(curveNames,curvesList):
        fwdRates = curve.fwd(fwdDates)
        plt.plot(years, fwdRates, label=name)
    plt.legend()
    plt.title("CC Fwd Rates")

    plt.figure()
    for name, curve in zip(curveNames,curvesList):
        fwdRates = curve.fwdRate(fwdDates, fwdDates2)
        plt.plot(years, fwdRates, label=name)
    plt.legend()
    plt.title("CC Fwd Rates")

    plt.figure()
    for name, curve in zip(curveNames,curvesList):
        dfs = curve.df(fwdDates)
        plt.plot(years, dfs, label=name)
    plt.legend()
    plt.title("Discount Factors")













###############################################################################


test_FinDiscountCurves()
testCases.compareTestCases()
