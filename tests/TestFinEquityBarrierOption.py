# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.models.FinProcessSimulator import FinProcessTypes
from financepy.models.FinProcessSimulator import FinGBMNumericalScheme
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierTypes
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierOption
from financepy.products.equity.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.finutils.FinDate import FinDate
import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinEquityBarrierOption():

    valueDate = FinDate(2015, 1, 1)
    expiryDate = FinDate(2016, 1, 1)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    optionType = FinEquityBarrierTypes.DOWN_AND_OUT_CALL

    drift = interestRate - dividendYield
    scheme = FinGBMNumericalScheme.NORMAL
    processType = FinProcessTypes.GBM
    discountCurve = FinFlatCurve(valueDate, interestRate)
    model = FinEquityModelBlackScholes(volatility)

    #######################################################################

    import time
    start = time.time()
    numObservationsPerYear = 100

    testCases.header(
        "Type",
        "K",
        "B",
        "S:",
        "Value:",
        "ValueMC",
        "Diff",
        "TIME")

    for optionType in FinEquityBarrierTypes:

        for stockPrice in range(80, 120, 10):

            B = 110.0
            K = 100.0

            option = FinEquityBarrierOption(
                expiryDate, K, optionType, B, numObservationsPerYear)
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model)
            start = time.time()
            modelParams = (stockPrice, drift, volatility, scheme)
            valueMC = option.valueMC(
                valueDate,
                stockPrice,
                discountCurve,
                processType,
                modelParams)

            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = valueMC - value

            testCases.print(
                optionType,
                K,
                B,
                stockPrice,
                value,
                valueMC,
                diff,
                timeElapsed)

    testCases.header(
        "Type",
        "K",
        "B",
        "S:",
        "Value:",
        "ValueMC",
        "Diff",
        "TIME")

    for stockPrice in range(80, 120, 10):

        B = 100.0
        K = 110.0

        option = FinEquityBarrierOption(
            expiryDate, K, optionType, B, numObservationsPerYear)
        value = option.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        start = time.time()
        modelParams = (stockPrice, drift, volatility, scheme)
        valueMC = option.valueMC(
            valueDate,
            stockPrice,
            discountCurve,
            processType,
            modelParams)
        end = time.time()
        timeElapsed = round(end - start, 3)
        diff = valueMC - value

        testCases.print(
            optionType,
            K,
            B,
            stockPrice,
            value,
            valueMC,
            diff,
            timeElapsed)

    end = time.time()


##########################################################################

    stockPrices = range(50, 150, 10)
    B = 105.0

    testCases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for optionType in FinEquityBarrierTypes:

        for stockPrice in stockPrices:

            barrierOption = FinEquityBarrierOption(
                expiryDate, 100.0, optionType, B, numObservationsPerYear)

            value = barrierOption.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model)
            delta = barrierOption.delta(
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model)
            vega = barrierOption.vega(
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model)
            theta = barrierOption.theta(
                valueDate,
                stockPrice,
                discountCurve,
                dividendYield,
                model)

            testCases.print(
                optionType,
                K,
                B,
                stockPrice,
                value,
                delta,
                vega,
                theta)


test_FinEquityBarrierOption()
testCases.compareTestCases()
