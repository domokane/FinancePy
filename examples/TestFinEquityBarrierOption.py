###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.models.process_simulator import FinProcessTypes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierTypes
from financepy.products.equity.FinEquityBarrierOption import FinEquityBarrierOption
from financepy.models.black_scholes import FinModelBlackScholes
from financepy.market.curves.FinDiscountCurveFlat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinEquityBarrierOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    optionType = FinEquityBarrierTypes.DOWN_AND_OUT_CALL

    drift = interestRate - dividendYield
    scheme = FinGBMNumericalScheme.NORMAL
    processType = FinProcessTypes.GBM

    discount_curve = DiscountCurveFlat(valuation_date, interestRate)
    dividendCurve = DiscountCurveFlat(valuation_date, dividendYield)

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
        for stock_price in range(80, 120, 10):

            B = 110.0
            K = 100.0

            option = FinEquityBarrierOption(
                expiry_date, K, optionType, B, numObservationsPerYear)
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)
            start = time.time()
            modelParams = (stock_price, drift, volatility, scheme)
            valueMC = option.valueMC(valuation_date,
                                     stock_price,
                                     discount_curve,
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
                stock_price,
                value,
                valueMC,
                diff,
                timeElapsed)

        for stock_price in range(80, 120, 10):

            B = 100.0
            K = 110.0

            option = FinEquityBarrierOption(
                expiry_date, K, optionType, B, numObservationsPerYear)
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)
            start = time.time()
            modelParams = (stock_price, drift, volatility, scheme)
            valueMC = option.valueMC(
                valuation_date,
                stock_price,
                discount_curve,
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
                stock_price,
                value,
                valueMC,
                diff,
                timeElapsed)

        end = time.time()

##########################################################################

    stock_prices = range(50, 150, 50)
    B = 105.0

    testCases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for optionType in FinEquityBarrierTypes:

        for stock_price in stock_prices:

            barrierOption = FinEquityBarrierOption(
                expiry_date, 100.0, optionType, B, numObservationsPerYear)

            value = barrierOption.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)
            delta = barrierOption.delta(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)
            vega = barrierOption.vega(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)
            theta = barrierOption.theta(
                valuation_date,
                stock_price,
                discount_curve,
                dividendCurve,
                model)

            testCases.print(
                optionType,
                K,
                B,
                stock_price,
                value,
                delta,
                vega,
                theta)

###############################################################################

test_FinEquityBarrierOption()
testCases.compareTestCases()
