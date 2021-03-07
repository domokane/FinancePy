###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

import numpy as np

from financepy.utils.global_types import FinOptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.error import FinError

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityVanillaOption():

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 7, 2015)
    stockPrice = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    numPathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for numPaths in numPathsList:

        callOption = EquityVanillaOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        start = time.time()
        value_mc = callOption.value_mc(valueDate, stockPrice, discountCurve,
                                     dividend_curve, model, numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, value_mc, duration)

###############################################################################

    stockPrices = range(80, 120, 10)
    numPaths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC", 
                     "CALL_VALUE_MC_SOBOL", "TIME")
    useSobol = True

    for stockPrice in stockPrices:

        callOption = EquityVanillaOption(expiryDate, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)

        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)

        start = time.time()

        useSobol = False
        value_mc1 = callOption.value_mc(valueDate, stockPrice, discountCurve,
                                      dividend_curve, model, numPaths, useSobol)

        useSobol = True
        value_mc2 = callOption.value_mc(valueDate, stockPrice, discountCurve,
                                      dividend_curve, model, numPaths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, value_mc1, value_mc2, duration)

###############################################################################

    stockPrices = range(80, 120, 10)
    numPaths = 100000

    testCases.header("NUMPATHS", "PUT_VALUE_BS", "PUT_VALUE_MC", 
                     "PUT_VALUE_MC_SOBOL", "TIME")

    for stockPrice in stockPrices:

        putOption = EquityVanillaOption(expiryDate, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valueDate, stockPrice, discountCurve,
                                dividend_curve, model)

        start = time.time()

        useSobol = False
        value_mc1 = putOption.value_mc(valueDate, stockPrice, discountCurve,
                                      dividend_curve, model, numPaths, useSobol)

        useSobol = True
        value_mc2 = putOption.value_mc(valueDate, stockPrice, discountCurve,
                                      dividend_curve, model, numPaths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, value_mc1, value_mc2, duration)

###############################################################################

    stockPrices = range(80, 120, 10)

    testCases.header("STOCK PRICE", "CALL_VALUE_BS", "CALL_DELTA_BS", 
                     "CALL_VEGA_BS", "CALL_THETA_BS", "CALL_RHO_BS")

    for stockPrice in stockPrices:

        callOption = EquityVanillaOption(expiryDate, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        delta = callOption.delta(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        vega = callOption.vega(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        theta = callOption.theta(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        rho = callOption.rho(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        testCases.print(stockPrice, value, delta, vega, theta, rho)

    ###########################################################################

    testCases.header("STOCK PRICE", "PUT_VALUE_BS", "PUT_DELTA_BS", 
                     "PUT_VEGA_BS", "PUT_THETA_BS", "PUT_RHO_BS")

    for stockPrice in stockPrices:
        
        putOption = EquityVanillaOption(expiryDate, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        delta = putOption.delta(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        vega = putOption.vega(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        theta = putOption.theta(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        rho = putOption.rho(valueDate, stockPrice, discountCurve,
                                 dividend_curve, model)
        testCases.print(stockPrice, value, delta, vega, theta, rho)


def testImpliedVolatility_NEW():


    valueDate = FinDate(1, 1, 2015)
    stockPrice = 100.0
    interestRate = 0.05
    dividendYield = 0.03
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    strikes = np.linspace(50, 150, 11)
    timesToExpiry = [0.003, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]    
    sigmas = np.arange(1, 100, 5) / 100.0
    option_types = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]

    testCases.header("OPT_TYPE", "TEXP", "STOCK_PRICE", "STRIKE", "INTRINSIC",
                     "VALUE", "INPUT_VOL", "IMPLIED_VOL")
    
    tol = 1e-5
    numTests = 0
    numFails = 0
    
    for vol in sigmas:

        model = BlackScholes(vol)

        for time_to_expiry in timesToExpiry:

            expiryDate = valueDate.addYears(time_to_expiry)

            for strike in strikes:

                for option_type in option_types:

                    option = EquityVanillaOption(expiryDate, strike, 
                                                    option_type)
                
                    value = option.value(valueDate, stockPrice, discountCurve, 
                                         dividend_curve, model)

                    intrinsic = option.intrinsic(valueDate, stockPrice,
                                             discountCurve, dividend_curve)

                    # I remove the cases where the time value is zero
                    # This is arbitrary but 1e-10 seems good enough to me
                    
                    impliedVol = -999

                    if value - intrinsic > 1e-10:

                        impliedVol = option.implied_volatility(valueDate,
                                                              stockPrice, 
                                                              discountCurve, 
                                                              dividend_curve,
                                                              value)
    
                    numTests += 1    
                        
                    errVol = np.abs(impliedVol - vol)
    
                    if errVol > tol:
    
                        testCases.print(option_type,
                                  time_to_expiry,
                                  stockPrice,
                                  strike, 
                                  intrinsic,
                                  value, 
                                  vol, 
                                  impliedVol)

                        # These fails include ones due to the zero time value    
                        numFails += 1
                            
                        testCases.print(option_type, time_to_expiry, stockPrice,
                                        strike,
                                        stockPrice, value, vol, impliedVol)

    assert numFails == 694, "Num Fails has changed."

#    print("Num Tests", numTests, "numFails", numFails)

###############################################################################

test_EquityVanillaOption()

start = time.time()
testImpliedVolatility_NEW()
end = time.time()
elapsed = end - start

#print("Elapsed:", elapsed)

testCases.compareTestCases()
