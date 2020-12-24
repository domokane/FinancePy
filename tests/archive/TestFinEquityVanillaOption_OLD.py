###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityVanillaOptionOLD import FinEquityVanillaOptionOLD
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityVanillaOptionFactored():

    valueDate = FinDate(2015, 1, 1)
    expiryDate = FinDate(2015, 7, 1)
    stockPrice = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    numPathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for numPaths in numPathsList:

        callOption = FinEquityVanillaOptionOLD(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        start = time.time()
        valueMC = callOption.valueMC(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)


###############################################################################

    stockPrices = range(80, 120, 10)
    numPaths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC", 
                     "CALL_VALUE_MC_SOBOL", "TIME")
    useSobol = True

    for stockPrice in stockPrices:

        callOption = FinEquityVanillaOptionOLD(expiryDate, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)

        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)

        start = time.time()

        useSobol = False
        valueMC1 = callOption.valueMC(valueDate, stockPrice, discountCurve,
                                      dividendYield, model, numPaths, useSobol)

        useSobol = True
        valueMC2 = callOption.valueMC(valueDate, stockPrice, discountCurve,
                                      dividendYield, model, numPaths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC1, valueMC2, duration)

###############################################################################

    stockPrices = range(80, 120, 10)
    numPaths = 100000

    testCases.header("NUMPATHS", "PUT_VALUE_BS", "PUT_VALUE_MC", 
                     "PUT_VALUE_MC_SOBOL", "TIME")

    for stockPrice in stockPrices:

        putOption = FinEquityVanillaOptionOLD(expiryDate, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valueDate, stockPrice, discountCurve,
                                dividendYield, model)

        start = time.time()

        useSobol = False
        valueMC1 = putOption.valueMC(valueDate, stockPrice, discountCurve,
                                      dividendYield, model, numPaths, useSobol)

        useSobol = True
        valueMC2 = putOption.valueMC(valueDate, stockPrice, discountCurve,
                                      dividendYield, model, numPaths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC1, valueMC2, duration)

###############################################################################

    stockPrices = range(80, 120, 10)

    testCases.header("STOCK PRICE", "CALL_VALUE_BS", "CALL_DELTA_BS", 
                     "CALL_VEGA_BS", "CALL_THETA_BS", "CALL_RHO_BS")

    for stockPrice in stockPrices:

        callOption = FinEquityVanillaOptionOLD(expiryDate, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        delta = callOption.delta(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        vega = callOption.vega(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        theta = callOption.theta(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        rho = callOption.rho(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        testCases.print(stockPrice, value, delta, vega, theta, rho)

    ###########################################################################

    testCases.header("STOCK PRICE", "PUT_VALUE_BS", "PUT_DELTA_BS", 
                     "PUT_VEGA_BS", "PUT_THETA_BS", "PUT_RHO_BS")

    for stockPrice in stockPrices:
        
        putOption = FinEquityVanillaOptionOLD(expiryDate, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        delta = putOption.delta(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        vega = putOption.vega(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        theta = putOption.theta(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        rho = putOption.rho(valueDate, stockPrice, discountCurve,
                                 dividendYield, model)
        testCases.print(stockPrice, value, delta, vega, theta, rho)

###############################################################################

    testCases.header("STOCK PRICE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    stockPrices = range(60, 150, 10)

    for stockPrice in stockPrices:
        callOption = FinEquityVanillaOptionOLD(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        impliedVol = callOption.impliedVolatility(
            valueDate, stockPrice, discountCurve, dividendYield, value)
        testCases.print(stockPrice, value, volatility, impliedVol)

###############################################################################


test_FinEquityVanillaOptionFactored()
testCases.compareTestCases()
