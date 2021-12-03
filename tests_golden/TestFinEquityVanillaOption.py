###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.error import FinError
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityVanillaOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    num_paths_list = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for num_paths in num_paths_list:

        call_option = EquityVanillaOption(
            expiry_date, 100.0, OptionTypes.EUROPEAN_CALL)

        value = call_option.value(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)
        start = time.time()

        value_mc = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                        dividend_curve, model, num_paths)
        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc, duration)


###############################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000

    testCases.header("NUMPATHS", "CALL_VALUE_BS", "CALL_VALUE_MC",
                     "CALL_VALUE_MC_SOBOL", "TIME")
    useSobol = True

    for stock_price in stock_prices:

        call_option = EquityVanillaOption(expiry_date, 100.0,
                                          OptionTypes.EUROPEAN_CALL)

        value = call_option.value(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)

        start = time.time()

        useSobol = False
        value_mc1 = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                         dividend_curve, model, num_paths, useSobol)

        useSobol = True
        value_mc2 = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                         dividend_curve, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc1, value_mc2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000

    testCases.header("NUMPATHS", "PUT_VALUE_BS", "PUT_VALUE_MC",
                     "PUT_VALUE_MC_SOBOL", "TIME")

    for stock_price in stock_prices:

        put_option = EquityVanillaOption(expiry_date, 100.0,
                                         OptionTypes.EUROPEAN_PUT)

        value = put_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_curve, model)

        start = time.time()

        useSobol = False
        value_mc1 = put_option.value_mc(valuation_date, stock_price, discount_curve,
                                        dividend_curve, model, num_paths, useSobol)

        useSobol = True
        value_mc2 = put_option.value_mc(valuation_date, stock_price, discount_curve,
                                        dividend_curve, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc1, value_mc2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)

    testCases.header("STOCK PRICE", "CALL_VALUE_BS", "CALL_DELTA_BS",
                     "CALL_VEGA_BS", "CALL_THETA_BS", "CALL_RHO_BS",
                     "CALL_VANNA_BS")

    for stock_price in stock_prices:

        call_option = EquityVanillaOption(expiry_date, 100.0,
                                          OptionTypes.EUROPEAN_CALL)
        value = call_option.value(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)
        delta = call_option.delta(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)
        vega = call_option.vega(valuation_date, stock_price, discount_curve,
                                dividend_curve, model)
        theta = call_option.theta(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)
        rho = call_option.rho(valuation_date, stock_price, discount_curve,
                              dividend_curve, model)
        vanna = call_option.vanna(valuation_date, stock_price, discount_curve,
                                  dividend_curve, model)
        testCases.print(stock_price, value, delta, vega, theta, rho, vanna)

    ###########################################################################

    testCases.header("STOCK PRICE", "PUT_VALUE_BS", "PUT_DELTA_BS",
                     "PUT_VEGA_BS", "PUT_THETA_BS", "PUT_RHO_BS", "PUT_VANNA_BS")

    for stock_price in stock_prices:

        put_option = EquityVanillaOption(expiry_date, 100.0,
                                         OptionTypes.EUROPEAN_PUT)

        value = put_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_curve, model)
        delta = put_option.delta(valuation_date, stock_price, discount_curve,
                                 dividend_curve, model)
        vega = put_option.vega(valuation_date, stock_price, discount_curve,
                               dividend_curve, model)
        theta = put_option.theta(valuation_date, stock_price, discount_curve,
                                 dividend_curve, model)
        rho = put_option.rho(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)
        vanna = put_option.vanna(valuation_date, stock_price, discount_curve,
                                 dividend_curve, model)
        testCases.print(stock_price, value, delta, vega, theta, rho, vanna)


def testImpliedVolatility_NEW():

    valuation_date = Date(1, 1, 2015)
    stock_price = 100.0
    interest_rate = 0.05
    dividend_yield = 0.03
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    strikes = np.linspace(50, 150, 11)
    timesToExpiry = [0.003, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]
    sigmas = np.arange(1, 100, 5) / 100.0
    option_types = [OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT]

    testCases.header("OPT_TYPE", "TEXP", "STOCK_PRICE", "STRIKE", "INTRINSIC",
                     "VALUE", "INPUT_VOL", "IMPLIED_VOL")

    tol = 1e-5
    numTests = 0
    numFails = 0

    for vol in sigmas:

        model = BlackScholes(vol)

        for time_to_expiry in timesToExpiry:

            expiry_date = valuation_date.add_years(time_to_expiry)

            for strike in strikes:

                for option_type in option_types:

                    option = EquityVanillaOption(expiry_date, strike,
                                                 option_type)

                    value = option.value(valuation_date, stock_price, discount_curve,
                                         dividend_curve, model)

                    intrinsic = option.intrinsic(valuation_date, stock_price,
                                                 discount_curve, dividend_curve)

                    # I remove the cases where the time value is zero
                    # This is arbitrary but 1e-10 seems good enough to me

                    impliedVol = -999

                    if value - intrinsic > 1e-10:

                        impliedVol = option.implied_volatility(valuation_date,
                                                               stock_price,
                                                               discount_curve,
                                                               dividend_curve,
                                                               value)

                    numTests += 1

                    errVol = np.abs(impliedVol - vol)

                    if errVol > tol:

                        testCases.print(option_type,
                                        time_to_expiry,
                                        stock_price,
                                        strike,
                                        intrinsic,
                                        value,
                                        vol,
                                        impliedVol)

                        # These fails include ones due to the zero time value
                        numFails += 1

                        testCases.print(option_type, time_to_expiry,
                                        stock_price,
                                        strike,
                                        stock_price, value, vol, impliedVol)

    assert numFails == 694, "Num Fails has changed."

#    print("Num Tests", numTests, "numFails", numFails)

###############################################################################


if 1 == 1:
    valuation_date = Date(30, 11, 2021)
    expiry_date = valuation_date.add_years(1)

    stock_price = 100
    volatility = 0.20
    model = BlackScholes(volatility)

    discount_curve = DiscountCurveFlat(valuation_date, 0.05)
    dividend_curve = DiscountCurveFlat(valuation_date, 0.0)

    call_option = EquityVanillaOption(
        expiry_date, 100.0, OptionTypes.EUROPEAN_CALL)

    value = call_option.value(expiry_date, 105.0, discount_curve,
                              dividend_curve, model)

    print(value)

else:
    test_EquityVanillaOption()

    start = time.time()
    testImpliedVolatility_NEW()
    end = time.time()
    elapsed = end - start

    # print("Elapsed:", elapsed)

    testCases.compareTestCases()
