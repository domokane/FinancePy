###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.fx.FinFXVanillaOption import FinFXVanillaOption
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################

def test_FinFXAmericanOption():

    # There is no FXAmericanOption class. It is embedded in the FXVanillaOption
    # class. This test just compares it to the European

    valueDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.030  # EUR
    ccy2CCRate = 0.025  # USD

    currencyPair = ccy1 + ccy2  # Always ccy1ccy2
    spotFXRate = 1.20
    strikeFXRate = 1.250
    volatility = 0.10

    domDiscountCurve = FinDiscountCurveFlat(valueDate, ccy2CCRate)
    forDiscountCurve = FinDiscountCurveFlat(valueDate, ccy1CCRate)

    model = FinModelBlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged

    testCases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spotFXRates = np.arange(50, 200, 10)/100.0

    for spotFXRate in spotFXRates:

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        currencyPair,
                                        FinOptionTypes.EUROPEAN_CALL,
                                        1000000,
                                        "USD")

        valueEuropean = callOption.value(valueDate,
                                         spotFXRate,
                                         domDiscountCurve,
                                         forDiscountCurve,
                                         model)['v']

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.AMERICAN_CALL,
                                        1000000,
                                        "USD")

        valueAmerican = callOption.value(valueDate,
                                         spotFXRate,
                                         domDiscountCurve,
                                         forDiscountCurve,
                                         model)['v']

        diff = (valueAmerican - valueEuropean)
        testCases.print(spotFXRate, valueEuropean, valueAmerican, diff)

    for spotFXRate in spotFXRates:

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_PUT,
                                        1000000,
                                        "USD")

        valueEuropean = callOption.value(valueDate,
                                         spotFXRate,
                                         domDiscountCurve,
                                         forDiscountCurve,
                                         model)['v']

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.AMERICAN_PUT,
                                        1000000,
                                        "USD")

        valueAmerican = callOption.value(valueDate,
                                         spotFXRate,
                                         domDiscountCurve,
                                         forDiscountCurve,
                                         model)['v']

        diff = (valueAmerican - valueEuropean)
        testCases.print(spotFXRate, valueEuropean, valueAmerican, diff)

###############################################################################

test_FinFXAmericanOption()
testCases.compareTestCases()
