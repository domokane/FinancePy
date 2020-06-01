# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.finutils.FinOptionTypes import FinOptionTypes
from financepy.products.fx.FinFXVanillaOption import FinFXVanillaOption
from financepy.products.fx.FinFXModelTypes import FinFXModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve

from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXOptionSABR():

    # UNFINISHED
    # There is no FXAmericanOption class. It is embedded in the FXVanillaOption
    # class. This test just compares it to the European

    valueDate = FinDate(13, 2, 2018)
    expiryDate = FinDate(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1CCRate = 0.030  # EUR
    ccy2CCRate = 0.025  # USD

    spotFXRate = 1.20
    strikeFXRate = 1.250
    volatility = 0.10

    notional = 1000000.0

    domDiscountCurve = FinFlatCurve(valueDate, ccy2CCRate)
    forDiscountCurve = FinFlatCurve(valueDate, ccy1CCRate)

    model = FinFXModelBlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged
    notional = 1000000.0

    testCases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spotFXRates = np.arange(50, 200, 10)/100.0

    for spotFXRate in spotFXRates:

        callOption = FinFXVanillaOption(expiryDate,
                                        strikeFXRate,
                                        "EURUSD",
                                        FinOptionTypes.EUROPEAN_CALL,
                                        notional,
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
        print("CALL %9.6f %9.6f %9.7f %10.8f" % (spotFXRate,
                                                 valueEuropean,
                                                 valueAmerican,
                                                 diff))

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
        print("PUT  %9.6f %9.6f %9.7f %10.8f" % (spotFXRate,
                                                 valueEuropean,
                                                 valueAmerican,
                                                 diff))

###############################################################################


test_FinFXOptionSABR()
testCases.compareTestCases()
