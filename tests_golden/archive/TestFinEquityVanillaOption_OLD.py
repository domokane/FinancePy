###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time

import sys
sys.path.append("..")

from financepy.utils.global_types import FinOptionTypes
from financepy.products.equity.EquityVanillaOptionOLD import EquityVanillaOptionOLD
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityVanillaOptionFactored():

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 7, 2015)
    stock_price = 100
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)

    num_paths_list = [10000, 20000, 40000, 80000, 160000, 320000]

    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    for num_paths in num_paths_list:

        call_option = EquityVanillaOptionOLD(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = call_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        start = time.time()
        value_mc = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                     dividend_yield, model, num_paths)
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

        call_option = EquityVanillaOptionOLD(expiry_date, 100.0,
                                            FinOptionTypes.EUROPEAN_CALL)

        value = call_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)

        start = time.time()

        useSobol = False
        value_mc1 = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                      dividend_yield, model, num_paths, useSobol)

        useSobol = True
        value_mc2 = call_option.value_mc(valuation_date, stock_price, discount_curve,
                                      dividend_yield, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc1, value_mc2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)
    num_paths = 100000

    testCases.header("NUMPATHS", "PUT_VALUE_BS", "PUT_VALUE_MC", 
                     "PUT_VALUE_MC_SOBOL", "TIME")

    for stock_price in stock_prices:

        put_option = EquityVanillaOptionOLD(expiry_date, 100.0,
                                           FinOptionTypes.EUROPEAN_PUT)

        value = put_option.value(valuation_date, stock_price, discount_curve,
                                dividend_yield, model)

        start = time.time()

        useSobol = False
        value_mc1 = put_option.value_mc(valuation_date, stock_price, discount_curve,
                                      dividend_yield, model, num_paths, useSobol)

        useSobol = True
        value_mc2 = put_option.value_mc(valuation_date, stock_price, discount_curve,
                                      dividend_yield, model, num_paths, useSobol)

        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc1, value_mc2, duration)

###############################################################################

    stock_prices = range(80, 120, 10)

    testCases.header("STOCK PRICE", "CALL_VALUE_BS", "CALL_DELTA_BS", 
                     "CALL_VEGA_BS", "CALL_THETA_BS", "CALL_RHO_BS")

    for stock_price in stock_prices:

        call_option = EquityVanillaOptionOLD(expiry_date, 100.0,
                                            FinOptionTypes.EUROPEAN_CALL)
        value = call_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        delta = call_option.delta(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        vega = call_option.vega(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        theta = call_option.theta(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        rho = call_option.rho(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        testCases.print(stock_price, value, delta, vega, theta, rho)

    ###########################################################################

    testCases.header("STOCK PRICE", "PUT_VALUE_BS", "PUT_DELTA_BS", 
                     "PUT_VEGA_BS", "PUT_THETA_BS", "PUT_RHO_BS")

    for stock_price in stock_prices:
        
        put_option = EquityVanillaOptionOLD(expiry_date, 100.0,
                                           FinOptionTypes.EUROPEAN_PUT)

        value = put_option.value(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        delta = put_option.delta(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        vega = put_option.vega(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        theta = put_option.theta(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        rho = put_option.rho(valuation_date, stock_price, discount_curve,
                                 dividend_yield, model)
        testCases.print(stock_price, value, delta, vega, theta, rho)

###############################################################################

    testCases.header("STOCK PRICE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    stock_prices = range(60, 150, 10)

    for stock_price in stock_prices:
        call_option = EquityVanillaOptionOLD(
            expiry_date, 100.0, FinOptionTypes.EUROPEAN_CALL)
        value = call_option.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_yield,
            model)
        impliedVol = call_option.implied_volatility(
            valuation_date, stock_price, discount_curve, dividend_yield, value)
        testCases.print(stock_price, value, volatility, impliedVol)

###############################################################################


test_EquityVanillaOptionFactored()
testCases.compareTestCases()
