###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

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

    testCases.header("STOCK PRICE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    stockPrices = range(60, 150, 10)

    for stockPrice in stockPrices:
        callOption = FinEquityVanillaOption(
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

    K = [10, 20, 50, 100, 150, 200, 250]
    T = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    sigma = [0.01, 0.10, 0.50, 1.0]
    optionTypes = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]
    stockPrice = 100
    HOURS_PER_YEAR = 365.25 * 24
    convergenceFailure = 0
    assertionFailure = 0
    noResult = 0
    successful = 0
    numberTests = 0
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

                    try:
                        impliedVol = option.impliedVolatility(
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
                            "Expected IV:", vol
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
                    finally:
                        numberTests += 1
    print("\nSuccessful:", successful)
    print("Convergence failure:", convergenceFailure)
    print("Inaccurate result:", assertionFailure)
    print("No result (price too low)", noResult)
    print("TOTAL:", numberTests)

###############################################################################

test_FinEquityVanillaOption()
testCases.compareTestCases()
