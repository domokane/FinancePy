###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

import numpy as np

from financepy.utils.global_types import FinOptionTypes
from financepy.products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import FinModelBlackScholes
from financepy.utils.date import Date
from financepy.utils.FinError import FinError

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityVanillaOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    model = FinModelBlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    num_pathsList = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for num_paths in num_pathsList:

        callOption = FinEquityVanillaOption(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        start = time.time()
        valueMC = callOption.valueMC(valuation_date, stock_price, discount_curve,
                                     dividendCurve, model, num_paths)
        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, valueMC, duration)

###############################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC", 
                     "CALL_VALUE_MC_SOBOL", "TIME")
    useSobol = True

    for stock_price in stock_prices:

        callOption = FinEquityVanillaOption(expiry_date, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)

        value = callOption.value(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)

        start = time.time()

        useSobol = False
        valueMC1 = callOption.valueMC(valuation_date, stock_price, discount_curve,
                                      dividendCurve, model, num_paths, useSobol)

        useSobol = True
        valueMC2 = callOption.valueMC(valuation_date, stock_price, discount_curve,
                                      dividendCurve, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, valueMC1, valueMC2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000

    testCases.header("NUMPATHS", "PUT_VALUE_BS", "PUT_VALUE_MC", 
                     "PUT_VALUE_MC_SOBOL", "TIME")

    for stock_price in stock_prices:

        putOption = FinEquityVanillaOption(expiry_date, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valuation_date, stock_price, discount_curve,
                                dividendCurve, model)

        start = time.time()

        useSobol = False
        valueMC1 = putOption.valueMC(valuation_date, stock_price, discount_curve,
                                      dividendCurve, model, num_paths, useSobol)

        useSobol = True
        valueMC2 = putOption.valueMC(valuation_date, stock_price, discount_curve,
                                      dividendCurve, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, valueMC1, valueMC2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)

    testCases.header("STOCK PRICE", "CALL_VALUE_BS", "CALL_DELTA_BS", 
                     "CALL_VEGA_BS", "CALL_THETA_BS", "CALL_RHO_BS")

    for stock_price in stock_prices:

        callOption = FinEquityVanillaOption(expiry_date, 100.0, 
                                            FinOptionTypes.EUROPEAN_CALL)
        value = callOption.value(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        delta = callOption.delta(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        vega = callOption.vega(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        theta = callOption.theta(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        rho = callOption.rho(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        testCases.print(stock_price, value, delta, vega, theta, rho)

    ###########################################################################

    testCases.header("STOCK PRICE", "PUT_VALUE_BS", "PUT_DELTA_BS", 
                     "PUT_VEGA_BS", "PUT_THETA_BS", "PUT_RHO_BS")

    for stock_price in stock_prices:
        
        putOption = FinEquityVanillaOption(expiry_date, 100.0, 
                                           FinOptionTypes.EUROPEAN_PUT)

        value = putOption.value(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        delta = putOption.delta(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        vega = putOption.vega(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        theta = putOption.theta(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        rho = putOption.rho(valuation_date, stock_price, discount_curve,
                                 dividendCurve, model)
        testCases.print(stock_price, value, delta, vega, theta, rho)


def testImpliedVolatility_NEW():


    valuation_date = Date(1, 1, 2015)
    stock_price = 100.0
    interestRate = 0.05
    dividendYield = 0.03
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    strikes = np.linspace(50, 150, 11)
    timesToExpiry = [0.003, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]    
    sigmas = np.arange(1, 100, 5) / 100.0
    optionTypes = [FinOptionTypes.EUROPEAN_CALL, FinOptionTypes.EUROPEAN_PUT]

    testCases.header("OPT_TYPE", "TEXP", "STOCK_PRICE", "STRIKE", "INTRINSIC",
                     "VALUE", "INPUT_VOL", "IMPLIED_VOL")
    
    tol = 1e-5
    numTests = 0
    numFails = 0
    
    for vol in sigmas:

        model = FinModelBlackScholes(vol)

        for timeToExpiry in timesToExpiry:     

            expiry_date = valuation_date.addYears(timeToExpiry)

            for strike in strikes:

                for optionType in optionTypes:

                    option = FinEquityVanillaOption(expiry_date, strike, 
                                                    optionType)
                
                    value = option.value(valuation_date, stock_price, discount_curve,
                                         dividendCurve, model)

                    intrinsic = option.intrinsic(valuation_date, stock_price,
                                             discount_curve, dividendCurve)

                    # I remove the cases where the time value is zero
                    # This is arbitrary but 1e-10 seems good enough to me
                    
                    impliedVol = -999

                    if value - intrinsic > 1e-10:

                        impliedVol = option.impliedVolatility(valuation_date, 
                                                              stock_price,
                                                              discount_curve, 
                                                              dividendCurve, 
                                                              value)
    
                    numTests += 1    
                        
                    errVol = np.abs(impliedVol - vol)
    
                    if errVol > tol:
    
                        testCases.print(optionType, 
                                  timeToExpiry, 
                                  stock_price,
                                  strike, 
                                  intrinsic,
                                  value, 
                                  vol, 
                                  impliedVol)

                        # These fails include ones due to the zero time value    
                        numFails += 1
                            
                        testCases.print(optionType, timeToExpiry, stock_price,
                                        strike,
                                        stock_price, value, vol, impliedVol)

    assert numFails == 694, "Num Fails has changed."

#    print("Num Tests", numTests, "numFails", numFails)

###############################################################################

test_FinEquityVanillaOption()

start = time.time()
testImpliedVolatility_NEW()
end = time.time()
elapsed = end - start

#print("Elapsed:", elapsed)

testCases.compareTestCases()
