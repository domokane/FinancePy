###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.models.FinProcessSimulator import FinProcessTypes
from financepy.models.FinProcessSimulator import FinGBMNumericalScheme
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierTypes
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierOption
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinEquityBarrierOption():

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    optionType = FinEquityBarrierTypes.DOWN_AND_OUT_CALL

    drift = interestRate - dividendYield
    scheme = FinGBMNumericalScheme.NORMAL
    processType = FinProcessTypes.GBM

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)    
    dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield)

    model = FinModelBlackScholes(volatility)

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
                dividendCurve,
                model)
            start = time.time()
            modelParams = (stockPrice, drift, volatility, scheme)
            valueMC = option.valueMC(valueDate,
                                     stockPrice,
                                     discountCurve,
                                     dividendCurve,
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

        for stockPrice in range(80, 120, 10):

            B = 100.0
            K = 110.0

            option = FinEquityBarrierOption(
                expiryDate, K, optionType, B, numObservationsPerYear)
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividendCurve,
                model)
            start = time.time()
            modelParams = (stockPrice, drift, volatility, scheme)
            valueMC = option.valueMC(
                valueDate,
                stockPrice,
                discountCurve,
                dividendCurve,
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

    stockPrices = range(50, 150, 50)
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
                dividendCurve,
                model)
            delta = barrierOption.delta(
                valueDate,
                stockPrice,
                discountCurve,
                dividendCurve,
                model)
            vega = barrierOption.vega(
                valueDate,
                stockPrice,
                discountCurve,
                dividendCurve,
                model)
            theta = barrierOption.theta(
                valueDate,
                stockPrice,
                discountCurve,
                dividendCurve,
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

###############################################################################

test_FinEquityBarrierOption()
testCases.compareTestCases()
