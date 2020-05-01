# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.equity.FinEquityCompoundOption import FinEquityCompoundOption
from financepy.products.equity.FinEquityOption import FinEquityOptionTypes
from financepy.products.equity.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.finutils.FinDate import FinDate

import sys
sys.path.append("..//..")


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

    model = FinEquityModelBlackScholes(volatility)
    discountCurve = FinFlatCurve(valueDate, interestRate)

    optionType1 = FinEquityOptionTypes.EUROPEAN_CALL
    optionType2 = FinEquityOptionTypes.EUROPEAN_PUT

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
        expiryDate2,
        k1,
        k2,
        optionType1,
        optionType2)
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
        values = cmpdOption.valueTree(
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
            FinEquityOptionTypes.EUROPEAN_CALL,
            FinEquityOptionTypes.EUROPEAN_PUT]:
        for optionType2 in [
                FinEquityOptionTypes.EUROPEAN_CALL,
                FinEquityOptionTypes.EUROPEAN_PUT]:

            cmpdOption = FinEquityCompoundOption(
                expiryDate1, expiryDate2, k1, k2, optionType1, optionType2)
            stockPrices = range(70, 100)

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

                values = cmpdOption.valueTree(valueDate, stockPrice,
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
