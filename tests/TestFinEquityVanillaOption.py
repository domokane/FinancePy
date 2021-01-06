###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

import numpy as np

from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinError import FinError

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityVanillaOption():

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

        callOption = FinEquityVanillaOption(
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

        callOption = FinEquityVanillaOption(expiryDate, 100.0, 
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

        putOption = FinEquityVanillaOption(expiryDate, 100.0, 
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

        callOption = FinEquityVanillaOption(expiryDate, 100.0, 
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
        
        putOption = FinEquityVanillaOption(expiryDate, 100.0, 
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

def testImpliedVolatility():


    valueDate = FinDate(1, 1, 2015)
    stockPrice = 100
    interestRate = 0.05
    dividendYield = 0.01
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    strikes = [10, 20, 50, 100, 150, 200]
    timesToExpiry = [0.003, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]    
    expiryDates = valueDate.addYears(timesToExpiry)
    sigmas = [0.01, 0.10, 0.50, 1.0]
    optionTypes = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]

    testCases.header("OPT_TYPE", "EXP_DATE", "STRIKE", "STOCK_PRICE",
                     "VALUE", "INPUT_VOL", "IMPLIED_VOL")
    
    tol = 1e-6
    numTests = 0
    numFails = 0
    
    for vol in sigmas:

        model = FinModelBlackScholes(vol)

        for expiryDate in expiryDates:     

            for strike in strikes:

                for optionType in optionTypes:

                    option = FinEquityVanillaOption(expiryDate, 100.0, 
                                                    optionType)
                
                    value = option.value(valueDate, stockPrice, discountCurve, 
                                         dividendYield, model)

                    impliedVol = option.impliedVolatility(valueDate, stockPrice, 
                                                      discountCurve, 
                                                      dividendYield, value)

                    numTests += 1    
                    if np.abs(impliedVol - vol) > tol:
                        numFails += 1

#                    print(optionType, expiryDate, strike, 
#                          stockPrice, value, vol, impliedVol)
            
                    testCases.print(optionType, expiryDate, strike, stockPrice, 
                                    value, vol, impliedVol)

    print("Num Tests", numTests, "numFails", numFails)

###############################################################################

    K = [10, 20, 50, 100, 150, 200]
    T = [0.003, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    sigma = [0.01, 0.10, 0.50, 1.0]
    optionTypes = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]
    stockPrice = 100
    HOURS_PER_YEAR = 365.25 * 24
    convergenceFailure = 0
    assertionFailure = 0
    noResult = 0
    successful = 0
    numberTests = 0
    totalElapsedTime = 0.
    for t in T:
        expDate = valueDate.addHours(t * HOURS_PER_YEAR)
        for k in K:
            for vol in sigma:
                bs_model = FinModelBlackScholes(vol)
                for type_ in optionTypes:
                    option = FinEquityVanillaOption(
                        expDate, k, type_)
                    value = option.value(
                        valueDate,
                        stockPrice,
                        discountCurve,
                        dividendYield,
                        bs_model)
                    if value < 1e-10:
                        continue
                    try:
                        start_time = time.time()
                        impliedVol = option.impliedVolatility_v2(
                            valueDate, stockPrice, discountCurve, dividendYield, value)
                        assert abs(impliedVol - vol) < 0.10
                    except FinError:
                        noResult += 1
                    except AssertionError:
                        assertionFailure += 1
                        print("-----------------")
                        print(f"Did not converge to expected value: {round(impliedVol, 2)} vs {vol}")
                        print(
                            "INPUTS\n",
                            "Type:", type_.name,
                            "ttm:", t,
                            "strike:", k,
                            "Expected IV:", vol,
                            "BS Price:", value
                        )
                    except RuntimeError:
                        import traceback
                        traceback.print_exc()
                        convergenceFailure += 1
                        print(
                            "INPUTS\n",
                            "Type:", type_.name,
                            "ttm:", t,
                            "strike:", k,
                            "Expected IV:", vol
                        )
                    except Exception:
                        import traceback
                        traceback.print_exc()
                        noResult += 1
                    else:
                        successful += 1
                        totalElapsedTime += time.time()-start_time
                    finally:
                        numberTests += 1

    print("\nSuccessful:", successful)
    print("Convergence failure:", convergenceFailure)
    print("Inaccurate result:", assertionFailure)
    print("No result (price too low)", noResult)
    print("TOTAL:", numberTests)
    print("Mean time:", 1e6 * totalElapsedTime/successful, "us")

###############################################################################

test_FinEquityVanillaOption()
start = time.time()
testImpliedVolatility()
end = time.time()
elapsed = end - start
# print("Elapsed:", elapsed)

testCases.compareTestCases()
