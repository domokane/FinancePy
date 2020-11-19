###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityDigitalOption import FinEquityDigitalOption, FinDigitalOptionTypes
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityDigitalOption():

    underlyingType = FinDigitalOptionTypes.CASH_OR_NOTHING

    valueDate = FinDate(2015, 1, 1)
    expiryDate = FinDate(2016, 1, 1)
    stockPrice = 100.0
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    model = FinModelBlackScholes(volatility)
    import time

    callOptionValues = []
    callOptionValuesMC = []
    numPathsList = [
        10000,
        20000,
        40000,
        80000,
        160000,
        320000,
        640000,
        1280000,
        2560000]

    testCases.header("NumLoops", "ValueBS", "ValueMC", "TIME")

    for numPaths in numPathsList:

        callOption = FinEquityDigitalOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, underlyingType)
        value = callOption.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        start = time.time()
        valueMC = callOption.valueMC(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model,
            numPaths)
        end = time.time()
        duration = end - start
        testCases.print(numPaths, value, valueMC, duration)

        callOptionValues.append(value)
        callOptionValuesMC.append(valueMC)

#    plt.figure(figsize=(10,8))
#    plt.plot(numPathsList, callOptionValues, color = 'b', label="Call Option")
#    plt.plot(numPathsList, callOptionValuesMC, color = 'r', label = "Call Option MC")
#    plt.xlabel("Num Loops")
#    plt.legend(loc='best')

##########################################################################

    stockPrices = range(50, 150, 50)
    callOptionValues = []
    callOptionDeltas = []
    callOptionVegas = []
    callOptionThetas = []

    for stockPrice in stockPrices:
        callOption = FinEquityDigitalOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_CALL, underlyingType)
        value = callOption.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        delta = callOption.delta(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vega = callOption.vega(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        theta = callOption.theta(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        callOptionValues.append(value)
        callOptionDeltas.append(delta)
        callOptionVegas.append(vega)
        callOptionThetas.append(theta)

    putOptionValues = []
    putOptionDeltas = []
    putOptionVegas = []
    putOptionThetas = []

    for stockPrice in stockPrices:
        putOption = FinEquityDigitalOption(
            expiryDate, 100.0, FinOptionTypes.EUROPEAN_PUT, underlyingType)
        value = putOption.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        delta = putOption.delta(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vega = putOption.vega(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        theta = putOption.theta(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        putOptionValues.append(value)
        putOptionDeltas.append(delta)
        putOptionVegas.append(vega)
        putOptionThetas.append(theta)

##########################################################################

test_FinEquityDigitalOption()
testCases.compareTestCases()
