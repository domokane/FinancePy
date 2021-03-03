###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import numpy as np
import matplotlib.pyplot as plt

from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

from financepy.models.FinModelVolatilityFns import FinVolFunctionTypes
from financepy.models.FinModelVolatilityFns import volFunctionClark

from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.models.FinModelOptionImpliedDbn import optionImpliedDbn

from financepy.market.volatility.FinFXVolSurface import FinFXVolSurface
from financepy.market.volatility.FinFXVolSurface import FinFXATMMethod
from financepy.market.volatility.FinFXVolSurface import FinFXDeltaMethod


from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinOptionImpliedDbn():

    if 1 == 1:

        # Example from Book extract by Iain Clark using Tables 3.3 and 3.4
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

#        fxMarket.checkCalibration(True)

        PLOT_GRAPHS = False
        if PLOT_GRAPHS:
            fxMarket.plotVolCurves()
 
        for iTenor in range(0, len(fxMarket._tenors)):
            
            F = fxMarket._F0T[iTenor]
            texp = fxMarket._texp[iTenor]

            startFX = F * 0.05
            endFX = F * 5.0

            numSteps = 10000
            dFX = (endFX - startFX)/ numSteps

            domDF = domDiscountCurve._df(texp)
            forDF = forDiscountCurve._df(texp)

            rd = -np.log(domDF) / texp
            rf = -np.log(forDF) / texp

            params = fxMarket._parameters[iTenor]

            strikes = []
            vols = []

            for iK in range(0, numSteps):
                strike = startFX + iK*dFX                
                vol = volFunctionClark(params, F, strike, texp)
                strikes.append(strike) 
                vols.append(vol)
            
            strikes = np.array(strikes)
            vols = np.array(vols)

            dbn = optionImpliedDbn(spotFXRate, texp, rd, rf, strikes, vols)
            
#            print("SUM:", dbn.sum())
#            plt.figure()
#            plt.plot(dbn._x, dbn._densitydx)

###############################################################################


test_FinOptionImpliedDbn()
testCases.compareTestCases()
