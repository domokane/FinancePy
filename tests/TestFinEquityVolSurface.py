###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import numpy as np

from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.volatility.FinEquityVolSurface import FinEquityVolSurface
from financepy.finutils.FinDate import FinDate
from financepy.models.FinModelVolatilityFns import FinVolFunctionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

import matplotlib.pyplot as plt

###############################################################################

PLOT_GRAPHS = False

###############################################################################
# TODO: ADD LOGGING TO TEST CASES
###############################################################################

def test_FinEquityVolSurfaceSVI(verboseCalibration):

    valueDate = FinDate(11, 1, 2021)

    stockPrice = 3800.0 # Check

    expiryDates = [FinDate(11, 2, 2021), FinDate(11, 3, 2021),
                   FinDate(11, 4, 2021), FinDate(11, 7, 2021), 
                   FinDate(11,10, 2021), FinDate(11, 1, 2022), 
                   FinDate(11, 1, 2023)]

    strikes = np.array([3037, 3418, 3608, 3703, 3798, 
                        3893, 3988, 4178, 4557])

    volSurface = [[42.94, 31.30, 25.88, 22.94, 19.72, 16.90, 15.31, 17.54, 25.67],
         [37.01, 28.25, 24.19, 21.93, 19.57, 17.45, 15.89, 15.34, 21.15],
         [34.68, 27.38, 23.82, 21.85, 19.83, 17.98, 16.52, 15.31, 18.94],
         [31.41, 26.25, 23.51, 22.05, 20.61, 19.25, 18.03, 16.01, 15.90],
         [29.91, 25.58, 23.21, 22.01, 20.83, 19.70, 18.62, 16.63, 14.94],
         [29.26, 25.24, 23.03, 21.91, 20.81, 19.73, 18.69, 16.76, 14.63],
         [27.59, 24.33, 22.72, 21.93, 21.17, 20.43, 19.71, 18.36, 16.26]]
        
    volSurface = np.array(volSurface)
    volSurface = volSurface / 100.0

    r = 0.020  # USD
    discountCurve = FinDiscountCurveFlat(valueDate, r)

    q = 0.010  # USD
    dividendCurve = FinDiscountCurveFlat(valueDate, q)

    volFunctionType = FinVolFunctionTypes.SVI

    equitySurface = FinEquityVolSurface(valueDate,
                                        stockPrice,
                                        discountCurve,
                                        dividendCurve,
                                        expiryDates,
                                        strikes,
                                        volSurface,
                                        volFunctionType)

#    tol = 1e-4
#    equitySurface.checkCalibration(False, tol)

    if 1==0: # PLOT_GRAPHS:

        equitySurface.plotVolCurves()

        plt.figure()

        mins = strikes[0] * 0.5
        maxs = strikes[-1] * 1.5

        dbns = equitySurface.impliedDbns(mins, maxs, 1000)

        for i in range(0, len(dbns)):
            expiryDateStr = str(equitySurface._expiryDates[i])
            plt.plot(dbns[i]._x, dbns[i]._densitydx, label = expiryDateStr)
            plt.title(volFunctionType)
            plt.legend()
            print("SUM:", dbns[i].sum())

    deltas = np.linspace(0.10, 0.90, 9)

    testCases.header("EXPIRY", "DELTA", "VOL", "STRIKE")
    for expiryDate in expiryDates:
        for delta in deltas:
            vol = equitySurface.volatilityFromDeltaDate(delta, expiryDate)
            testCases.print(expiryDate, delta, vol[0], vol[1])

###############################################################################

def test_FinEquityVolSurfaceSSVI(verboseCalibration):

    valueDate = FinDate(11, 1, 2021)

    stockPrice = 3800.0 # Check

    expiryDates = [FinDate(11, 2, 2021), FinDate(11, 3, 2021),
                   FinDate(11, 4, 2021), FinDate(11, 7, 2021), 
                   FinDate(11,10, 2021), FinDate(11, 1, 2022), 
                   FinDate(11, 1, 2023)]

    strikes = np.array([3037, 3418, 3608, 3703, 3798, 
                        3893, 3988, 4178, 4557])

    volSurface = [[42.94, 31.30, 25.88, 22.94, 19.72, 16.90, 15.31, 17.54, 25.67],
         [37.01, 28.25, 24.19, 21.93, 19.57, 17.45, 15.89, 15.34, 21.15],
         [34.68, 27.38, 23.82, 21.85, 19.83, 17.98, 16.52, 15.31, 18.94],
         [31.41, 26.25, 23.51, 22.05, 20.61, 19.25, 18.03, 16.01, 15.90],
         [29.91, 25.58, 23.21, 22.01, 20.83, 19.70, 18.62, 16.63, 14.94],
         [29.26, 25.24, 23.03, 21.91, 20.81, 19.73, 18.69, 16.76, 14.63],
         [27.59, 24.33, 22.72, 21.93, 21.17, 20.43, 19.71, 18.36, 16.26]]
        
    volSurface = np.array(volSurface)
    volSurface = volSurface / 100.0

    r = 0.020  # USD
    discountCurve = FinDiscountCurveFlat(valueDate, r)

    q = 0.010  # USD
    dividendCurve = FinDiscountCurveFlat(valueDate, q)

    volFunctionType = FinVolFunctionTypes.SSVI

    equitySurface = FinEquityVolSurface(valueDate,
                                        stockPrice,
                                        discountCurve,
                                        dividendCurve,
                                        expiryDates,
                                        strikes,
                                        volSurface,
                                        volFunctionType)

#    tol = 1e-4
#    equitySurface.checkCalibration(False, tol)

    if 1==0: # PLOT_GRAPHS:

        equitySurface.plotVolCurves()

        plt.figure()

        mins = strikes[0] * 0.5
        maxs = strikes[-1] * 1.5

        dbns = equitySurface.impliedDbns(mins, maxs, 1000)

        for i in range(0, len(dbns)):
            expiryDateStr = str(equitySurface._expiryDates[i])
            plt.plot(dbns[i]._x, dbns[i]._densitydx, label = expiryDateStr)
            plt.title(volFunctionType)
            plt.legend()
            print("SUM:", dbns[i].sum())

    deltas = np.linspace(0.10, 0.90, 9)

    testCases.header("EXPIRY", "DELTA", "VOL", "STRIKE")
    for expiryDate in expiryDates:
        for delta in deltas:
            vol = equitySurface.volatilityFromDeltaDate(delta, expiryDate)
            testCases.print(expiryDate, delta, vol[0], vol[1])


###############################################################################

import time

if __name__ == '__main__':

    start = time.time()

    verboseCalibration = False

#    test_FinEquityVolSurfaceSVI(verboseCalibration)
    test_FinEquityVolSurfaceSSVI(verboseCalibration)
    
    end = time.time()
    
    elapsed = end - start
    print("Elapsed Time:", elapsed)
    testCases.compareTestCases()
