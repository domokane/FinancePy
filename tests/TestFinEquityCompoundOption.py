###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import FinOptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityCompoundOption():

    valueDate = FinDate(1, 1, 2015)
    expiryDate1 = FinDate(1, 1, 2017)
    expiryDate2 = FinDate(1, 1, 2018)
    k1 = 5.0
    k2 = 95.0
    stockPrice = 85.0
    volatility = 0.15
    interestRate = 0.035
    dividendYield = 0.01

    model = BlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    numStepsList = [100, 200, 500, 1000, 2000, 5000]

    ###########################################################################

    stockPrice = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "TreeSteps", "Exact", "TreeValue")

    for option_type1 in [
            FinOptionTypes.EUROPEAN_CALL,
            FinOptionTypes.EUROPEAN_PUT]:
        for option_type2 in [
                FinOptionTypes.EUROPEAN_CALL,
                FinOptionTypes.EUROPEAN_PUT]:

            cmpdOption = EquityCompoundOption(expiryDate1, option_type1, k1,
                                                 expiryDate2, option_type2, k2)

            for numSteps in numStepsList:
        
                value = cmpdOption.value(valueDate, stockPrice, discountCurve,
                                         dividend_curve, model)

                values = cmpdOption._valueTree(valueDate, stockPrice, discountCurve,
                                               dividend_curve, model, numSteps)
        
                testCases.print(option_type1, option_type2, k1, k2, stockPrice,
                                numSteps, value, values[0])

    ###########################################################################

    stockPrice = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "TreeSteps", "Exact", "TreeValue")

    for option_type1 in [
            FinOptionTypes.AMERICAN_CALL,
            FinOptionTypes.AMERICAN_PUT]:
        for option_type2 in [
                FinOptionTypes.AMERICAN_CALL,
                FinOptionTypes.AMERICAN_PUT]:

            cmpdOption = EquityCompoundOption(expiryDate1, option_type1, k1,
                                                 expiryDate2, option_type2, k2)

            for numSteps in numStepsList:
        
                value = cmpdOption.value(valueDate, stockPrice, discountCurve,
                                         dividend_curve, model, numSteps)

                values = cmpdOption._valueTree(valueDate, stockPrice, discountCurve,
                                               dividend_curve, model, numSteps)
        
                testCases.print(option_type1, option_type2, k1, k2, stockPrice,
                                numSteps, value, values[0])

    ###########################################################################

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "Exact", "TreeSteps",
                     "TreeValue", "Diff", "DELTA", "GAMMA", "THETA")

    for option_type1 in [
            FinOptionTypes.EUROPEAN_CALL,
            FinOptionTypes.EUROPEAN_PUT]:
        for option_type2 in [
                FinOptionTypes.EUROPEAN_CALL,
                FinOptionTypes.EUROPEAN_PUT]:

            cmpdOption = EquityCompoundOption(
                expiryDate1, option_type1, k1,
                expiryDate2, option_type2, k2)
            stockPrices = range(70, 100, 10)

            for stockPrice in stockPrices:
                value = cmpdOption.value(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividend_curve,
                    model)
                delta = cmpdOption.delta(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividend_curve,
                    model)
                vega = cmpdOption.vega(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividend_curve,
                    model)
                theta = cmpdOption.theta(
                    valueDate,
                    stockPrice,
                    discountCurve,
                    dividend_curve,
                    model)

                values = cmpdOption._valueTree(valueDate, stockPrice,
                                               discountCurve, dividend_curve,
                                               model)

                diff = value - values[0]

                testCases.print(
                    option_type1,
                    option_type2,
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


test_EquityCompoundOption()
testCases.compareTestCases()
