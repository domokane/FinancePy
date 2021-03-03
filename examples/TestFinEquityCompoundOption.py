###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityCompoundOption import FinEquityCompoundOption
from financepy.utils.FinGlobalTypes import FinOptionTypes
from financepy.models.black_scholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinEquityCompoundOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date1 = Date(1, 1, 2017)
    expiry_date2 = Date(1, 1, 2018)
    k1 = 5.0
    k2 = 95.0
    stock_price = 85.0
    volatility = 0.15
    interestRate = 0.035
    dividendYield = 0.01

    model = FinModelBlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

    num_stepsList = [100, 200, 500, 1000, 2000, 5000]

    ###########################################################################

    stock_price = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "TreeSteps", "Exact", "TreeValue")

    for optionType1 in [
            FinOptionTypes.EUROPEAN_CALL,
            FinOptionTypes.EUROPEAN_PUT]:
        for optionType2 in [
                FinOptionTypes.EUROPEAN_CALL,
                FinOptionTypes.EUROPEAN_PUT]:

            cmpdOption = FinEquityCompoundOption(expiry_date1, optionType1, k1,
                                                 expiry_date2, optionType2, k2)

            for num_steps in num_stepsList:
        
                value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                                         dividendCurve, model)

                values = cmpdOption._valueTree(valuation_date, stock_price, discount_curve,
                                               dividendCurve, model, num_steps)
        
                testCases.print(optionType1, optionType2, k1, k2, stock_price,
                                num_steps, value, values[0])

    ###########################################################################

    stock_price = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "TreeSteps", "Exact", "TreeValue")

    for optionType1 in [
            FinOptionTypes.AMERICAN_CALL,
            FinOptionTypes.AMERICAN_PUT]:
        for optionType2 in [
                FinOptionTypes.AMERICAN_CALL,
                FinOptionTypes.AMERICAN_PUT]:

            cmpdOption = FinEquityCompoundOption(expiry_date1, optionType1, k1,
                                                 expiry_date2, optionType2, k2)

            for num_steps in num_stepsList:
        
                value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                                         dividendCurve, model, num_steps)

                values = cmpdOption._valueTree(valuation_date, stock_price, discount_curve,
                                               dividendCurve, model, num_steps)
        
                testCases.print(optionType1, optionType2, k1, k2, stock_price,
                                num_steps, value, values[0])

    ###########################################################################

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "Exact", "TreeSteps",
                     "TreeValue", "Diff", "DELTA", "GAMMA", "THETA")

    for optionType1 in [
            FinOptionTypes.EUROPEAN_CALL,
            FinOptionTypes.EUROPEAN_PUT]:
        for optionType2 in [
                FinOptionTypes.EUROPEAN_CALL,
                FinOptionTypes.EUROPEAN_PUT]:

            cmpdOption = FinEquityCompoundOption(
                expiry_date1, optionType1, k1,
                expiry_date2, optionType2, k2)
            stock_prices = range(70, 100, 10)

            for stock_price in stock_prices:
                value = cmpdOption.value(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividendCurve,
                    model)
                delta = cmpdOption.delta(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividendCurve,
                    model)
                vega = cmpdOption.vega(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividendCurve,
                    model)
                theta = cmpdOption.theta(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividendCurve,
                    model)

                values = cmpdOption._valueTree(valuation_date, stock_price,
                                               discount_curve, dividendCurve,
                                               model)

                diff = value - values[0]

                testCases.print(
                    optionType1,
                    optionType2,
                    k1,
                    k2,
                    stock_price,
                    value,
                    num_steps,
                    values[0],
                    diff,
                    delta,
                    vega,
                    theta)

##########################################################################


test_FinEquityCompoundOption()
testCases.compareTestCases()
