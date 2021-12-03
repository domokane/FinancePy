###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.products.equity.equity_barrier_option import EquityBarrierTypes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityBarrierOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    option_type = EquityBarrierTypes.DOWN_AND_OUT_CALL

    drift = interest_rate - dividend_yield
    scheme = FinGBMNumericalScheme.NORMAL
    process_type = ProcessTypes.GBM

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

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
        for stock_price in range(80, 120, 10):

            B = 110.0
            K = 100.0

            option = EquityBarrierOption(
                expiry_date, K, option_type, B, num_observations_per_year)
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)
            start = time.time()
            model_params = (stock_price, drift, volatility, scheme)
            value_mc = option.value_mc(valuation_date,
                                       stock_price,
                                       discount_curve,
                                       dividend_curve,
                                       process_type,
                                       model_params)

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(
                option_type,
                K,
                B,
                stock_price,
                value,
                value_mc,
                diff,
                time_elapsed)

        for stock_price in range(80, 120, 10):

            B = 100.0
            K = 110.0

            option = EquityBarrierOption(
                expiry_date, K, option_type, B, num_observations_per_year)
            value = option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)
            start = time.time()
            model_params = (stock_price, drift, volatility, scheme)
            value_mc = option.value_mc(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                process_type,
                model_params)
            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            testCases.print(
                option_type,
                K,
                B,
                stock_price,
                value,
                value_mc,
                diff,
                time_elapsed)

        end = time.time()

##########################################################################

    stock_prices = range(50, 150, 50)
    B = 105.0

    testCases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for option_type in EquityBarrierTypes:

        for stock_price in stock_prices:

            barrier_option = EquityBarrierOption(
                expiry_date, 100.0, option_type, B, num_observations_per_year)

            value = barrier_option.value(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)
            delta = barrier_option.delta(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)
            vega = barrier_option.vega(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)
            theta = barrier_option.theta(
                valuation_date,
                stock_price,
                discount_curve,
                dividend_curve,
                model)

            testCases.print(
                option_type,
                K,
                B,
                stock_price,
                value,
                delta,
                vega,
                theta)

###############################################################################


test_EquityBarrierOption()
testCases.compareTestCases()
