###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.models.FinProcessSimulator import FinProcessTypes
from financepy.models.FinProcessSimulator import FinGBMNumericalScheme
from financepy.products.equity.equity_barrier_option import EquityBarrierTypes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_EquityBarrierOption():

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 1, 2016)
    stockPrice = 100.0
    volatility = 0.20
    interestRate = 0.05
    dividendYield = 0.02
    option_type = EquityBarrierTypes.DOWN_AND_OUT_CALL

    drift = interestRate - dividendYield
    scheme = FinGBMNumericalScheme.NORMAL
    process_type = FinProcessTypes.GBM

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)    
    dividend_curve = FinDiscountCurveFlat(valueDate, dividendYield)

    model = BlackScholes(volatility)

    #######################################################################

    import time
    start = time.time()
    num_observations_per_year = 100

    testCases.header(
        "Type",
        "K",
        "B",
        "S:",
        "Value:",
        "ValueMC",
        "Diff",
        "TIME")

    for option_type in EquityBarrierTypes:
        for stockPrice in range(80, 120, 10):

            B = 110.0
            K = 100.0

            option = EquityBarrierOption(
                expiryDate, K, option_type, B, num_observations_per_year)
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)
            start = time.time()
            model_params = (stockPrice, drift, volatility, scheme)
            value_mc = option.value_mc(valueDate,
                                     stockPrice,
                                     discountCurve,
                                     dividend_curve,
                                     process_type,
                                     model_params)

            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(
                option_type,
                K,
                B,
                stockPrice,
                value,
                value_mc,
                diff,
                timeElapsed)

        for stockPrice in range(80, 120, 10):

            B = 100.0
            K = 110.0

            option = EquityBarrierOption(
                expiryDate, K, option_type, B, num_observations_per_year)
            value = option.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)
            start = time.time()
            model_params = (stockPrice, drift, volatility, scheme)
            value_mc = option.value_mc(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                process_type,
                model_params)
            end = time.time()
            timeElapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(
                option_type,
                K,
                B,
                stockPrice,
                value,
                value_mc,
                diff,
                timeElapsed)

        end = time.time()

##########################################################################

    stockPrices = range(50, 150, 50)
    B = 105.0

    testCases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for option_type in EquityBarrierTypes:

        for stockPrice in stockPrices:

            barrierOption = EquityBarrierOption(
                expiryDate, 100.0, option_type, B, num_observations_per_year)

            value = barrierOption.value(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)
            delta = barrierOption.delta(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)
            vega = barrierOption.vega(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)
            theta = barrierOption.theta(
                valueDate,
                stockPrice,
                discountCurve,
                dividend_curve,
                model)

            testCases.print(
                option_type,
                K,
                B,
                stockPrice,
                value,
                delta,
                vega,
                theta)

###############################################################################

test_EquityBarrierOption()
testCases.compareTestCases()
