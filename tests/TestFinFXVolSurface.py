###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.market.volatility.FinFXVolSurface import FinFXVolSurface
from financepy.market.volatility.FinFXVolSurface import FinFXATMMethod
from financepy.market.volatility.FinFXVolSurface import FinFXDeltaMethod
from financepy.finutils.FinDate import FinDate
from financepy.models.FinModelVolatilityFns import FinVolFunctionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

import matplotlib.pyplot as plt

###############################################################################

PLOT_GRAPHS = False

###############################################################################


def test_FinFXMktVolSurface1(verboseCalibration):

    ###########################################################################

    if 1 == 1:

        # Example from Book extract by Iain Clarke using Tables 3.3 and 3.4
        # print("EURUSD EXAMPLE CLARK")

        valueDate = FinDate(10, 4, 2020)

        forName = "EUR"
        domName = "USD"
        forCCRate = 0.03460  # EUR
        domCCRate = 0.02940  # USD

        domDiscountCurve = FinDiscountCurveFlat(valueDate, domCCRate)
        forDiscountCurve = FinDiscountCurveFlat(valueDate, forCCRate)

        currencyPair = forName + domName
        spotFXRate = 1.3465

        tenors = ['1M', '2M', '3M', '6M', '1Y', '2Y']
        atmVols = [21.00, 21.00, 20.750, 19.400, 18.250, 17.677]
        marketStrangle25DeltaVols = [0.65, 0.75, 0.85, 0.90, 0.95, 0.85]
        riskReversal25DeltaVols = [-0.20, -0.25, -0.30, -0.50, -0.60, -0.562]

        notionalCurrency = forName

        atmMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL
        deltaMethod = FinFXDeltaMethod.SPOT_DELTA
        volFunctionType = FinVolFunctionTypes.CLARK

        fxMarket = FinFXVolSurface(valueDate,
                                   spotFXRate,
                                   currencyPair,
                                   notionalCurrency,
                                   domDiscountCurve,
                                   forDiscountCurve,
                                   tenors,
                                   atmVols,
                                   marketStrangle25DeltaVols,
                                   riskReversal25DeltaVols,
                                   atmMethod,
                                   deltaMethod, 
                                   volFunctionType)

        fxMarket.checkCalibration(verboseCalibration)

        # EXPLORE AND TEST DIFFERENT CATEGORICAL PARAMETERS
        # for atmMethod in FinFXATMMethod:
        #     for deltaMethod in FinFXDeltaMethod:
        #         for volFunctionType in FinVolFunctionTypes:
    
        #             fxMarket = FinFXVolSurface(valueDate,
        #                                        spotFXRate,
        #                                        currencyPair,
        #                                        notionalCurrency,
        #                                        domDiscountCurve,
        #                                        forDiscountCurve,
        #                                        tenors,
        #                                        atmVols,
        #                                        marketStrangle25DeltaVols,
        #                                        riskReversal25DeltaVols,
        #                                        atmMethod,
        #                                        deltaMethod, 
        #                                        volFunctionType)
        
        #             fxMarket.checkCalibration(verboseCalibration)

        if PLOT_GRAPHS:

            fxMarket.plotVolCurves()

            dbns = fxMarket.impliedDbns(0.00001, 5.0, 10000)

            for i in range(0, len(dbns)):
                plt.plot(dbns[i]._x, dbns[i]._densitydx)
                plt.title(volFunctionType)
                print("SUM:", dbns[i].sum())

    ###########################################################################

def test_FinFXMktVolSurface2(verboseCalibration):

        #print("==============================================================")

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
        # print("EURJPY EXAMPLE CLARK")

        valueDate = FinDate(10, 4, 2020)

        forName = "EUR"
        domName = "JPY"
        forCCRate = 0.0294  # EUR
        domCCRate = 0.0171  # USD

        domDiscountCurve = FinDiscountCurveFlat(valueDate, domCCRate)
        forDiscountCurve = FinDiscountCurveFlat(valueDate, forCCRate)

        currencyPair = forName + domName
        spotFXRate = 90.72

        tenors = ['1M', '2M', '3M', '6M', '1Y', '2Y']
        atmVols = [21.50, 20.50, 19.85, 18.00, 15.95, 14.009]
        marketStrangle25DeltaVols = [0.35, 0.325, 0.300, 0.225, 0.175, 0.100]
        riskReversal25DeltaVols = [-8.350, -8.650, -8.950, -9.250, -9.550, -9.500]

        notionalCurrency = forName

        atmMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL_PREM_ADJ
        deltaMethod = FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ

        fxMarket = FinFXVolSurface(valueDate,
                                   spotFXRate,
                                   currencyPair,
                                   notionalCurrency,
                                   domDiscountCurve,
                                   forDiscountCurve,
                                   tenors,
                                   atmVols,
                                   marketStrangle25DeltaVols,
                                   riskReversal25DeltaVols,
                                   atmMethod,
                                   deltaMethod)

        fxMarket.checkCalibration(verboseCalibration)

        if PLOT_GRAPHS:
            fxMarket.plotVolCurves()

#   print("==================================================================")

#   ###########################################################################

def test_FinFXMktVolSurface3(verboseCalibration):

        # EURUSD Example from Paper by Uwe Wystup using Tables 4
#        print("EURUSD EXAMPLE WYSTUP")

        valueDate = FinDate(20, 1, 2009)

        forName = "EUR"
        domName = "USD"
        forCCRate = 0.020113  # EUR
        domCCRate = 0.003525  # USD

        domDiscountCurve = FinDiscountCurveFlat(valueDate, domCCRate)
        forDiscountCurve = FinDiscountCurveFlat(valueDate, forCCRate)

        currencyPair = forName + domName
        spotFXRate = 1.3088

        tenors = ['1M']
        atmVols = [21.6215]
        marketStrangle25DeltaVols = [0.7375]
        riskReversal25DeltaVols = [-0.50]

        notionalCurrency = forName

        atmMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL
        deltaMethod = FinFXDeltaMethod.SPOT_DELTA

        fxMarket = FinFXVolSurface(valueDate,
                                   spotFXRate,
                                   currencyPair,
                                   notionalCurrency,
                                   domDiscountCurve,
                                   forDiscountCurve,
                                   tenors,
                                   atmVols,
                                   marketStrangle25DeltaVols,
                                   riskReversal25DeltaVols,
                                   atmMethod,
                                   deltaMethod)

        fxMarket.checkCalibration(verboseCalibration)

        if PLOT_GRAPHS:
            fxMarket.plotVolCurves()

    ###########################################################################

def test_FinFXMktVolSurface4(verboseCalibration):

        # USDJPY Example from Paper by Uwe Wystup using Tables 4
#        print("USDJPY EXAMPLE WYSTUP")

        valueDate = FinDate(20, 1, 2009)

        forName = "USD"
        domName = "JPY"
        forCCRate = 0.003525  # USD
        domCCRate = 0.0042875  # JPY

        domDiscountCurve = FinDiscountCurveFlat(valueDate, domCCRate)
        forDiscountCurve = FinDiscountCurveFlat(valueDate, forCCRate)

        currencyPair = forName + domName
        spotFXRate = 90.68

        tenors = ['1M']
        atmVols = [21.00]
        marketStrangle25DeltaVols = [0.184]
        riskReversal25DeltaVols = [-5.30]

        notionalCurrency = forName

        atmMethod = FinFXATMMethod.FWD_DELTA_NEUTRAL
        deltaMethod = FinFXDeltaMethod.SPOT_DELTA_PREM_ADJ

        fxMarket = FinFXVolSurface(valueDate,
                                   spotFXRate,
                                   currencyPair,
                                   notionalCurrency,
                                   domDiscountCurve,
                                   forDiscountCurve,
                                   tenors,
                                   atmVols,
                                   marketStrangle25DeltaVols,
                                   riskReversal25DeltaVols,
                                   atmMethod,
                                   deltaMethod)

        fxMarket.checkCalibration(verboseCalibration)

        if PLOT_GRAPHS:
            fxMarket.plotVolCurves()

    #    testCases.header("value", "delta")
    #    testCases.print(value, delta)

###############################################################################

import time

if __name__ == '__main__':

    start = time.time()

    verboseCalibration = False

    test_FinFXMktVolSurface1(verboseCalibration)
    test_FinFXMktVolSurface2(verboseCalibration)
    test_FinFXMktVolSurface3(verboseCalibration)
    test_FinFXMktVolSurface4(verboseCalibration)
    
    end = time.time()
    
    elapsed = end - start
#    print("Elapsed Time:", elapsed)
    testCases.compareTestCases()
