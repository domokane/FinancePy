###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityCompoundOption import FinEquityCompoundOption
from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityCompoundOption():

    valueDate = FinDate(2015, 1, 1)
    expiryDate1 = FinDate(2017, 1, 1)
    expiryDate2 = FinDate(2018, 1, 1)
    k1 = 5.0
    k2 = 95.0
    stockPrice = 85.0
    volatility = 0.15
    interestRate = 0.035
    dividendYield = 0.01

    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    optionType1 = FinOptionTypes.EUROPEAN_CALL
    optionType2 = FinOptionTypes.EUROPEAN_PUT

    numStepsList = [200,
                    300,
                    400,
                    500,
                    600,
                    700,
                    800,
                    900,
                    1000,
                    2000,
                    3000,
                    4000,
                    5000]

    cmpdOption = FinEquityCompoundOption(
        expiryDate1,
        optionType1,
        k1,
        expiryDate2,
        optionType2,
        k2)
    stockPrice = 85.0

    testCases.header(
        "TYPE1",
        "TYPE2",
        "K1",
        "K2",
        "S",
        "Exact",
        "TreeSteps",
        "TreeValue")
    for numSteps in numStepsList:

        value = cmpdOption.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        values = cmpdOption._valueTree(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model,
            numSteps)
        testCases.print(
            optionType1,
            optionType2,
            k1,
            k2,
            stockPrice,
            value,
            numSteps,
            values[0])

    testCases.header(
        "TYPE1",
        "TYPE2",
        "K1",
        "K2",
        "S",
        "Exact",
        "TreeSteps",
        "TreeValue",
        "Diff",
        "DELTA",
        "GAMMA",
        "THETA")

    for optionType1 in [
            FinOptionTypes.EUROPEAN_CALL,
            FinOptionTypes.EUROPEAN_PUT]:
        for optionType2 in [
                FinOptionTypes.EUROPEAN_CALL,
                FinOptionTypes.EUROPEAN_PUT]:

            cmpdOption = FinEquityCompoundOption(
                expiryDate1, optionType1, k1, expiryDate2, optionType2, k2)
            stockPrices = range(70, 100, 5)

            for stockPrice in stockPrices:
                value = cmpdOption.value(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividendYield,
                    model)
                delta = cmpdOption.delta(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividendYield,
                    model)
                vega = cmpdOption.vega(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividendYield,
                    model)
                theta = cmpdOption.theta(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividendYield,
                    model)

                values = cmpdOption._valueTree(valueDate, stockPrice,
                                               discountCurve, dividendYield,
                                               model)

                diff = value - values[0]

                testCases.print(
                    optionType1,
                    optionType2,
                    k1,
                    k2,
                    stockPrice,
                    value,
                    numSteps,
                    values[0],
                    diff,
                    delta,
                    vega,
                    theta)

##########################################################################

test_FinEquityCompoundOption()
testCases.compareTestCases()
