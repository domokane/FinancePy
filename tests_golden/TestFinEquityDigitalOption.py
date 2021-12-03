###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_digital_option import EquityDigitalOption, FinDigitalOptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityDigitalOption():

    underlying_type = FinDigitalOptionTypes.CASH_OR_NOTHING

    valuation_date = Date(1, 1, 2015)
    expiry_date = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.30
    interest_rate = 0.05
    dividend_yield = 0.01
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    model = BlackScholes(volatility)
    import time

    call_option_values = []
    call_option_valuesMC = []
    num_paths_list = [
        10000,
        20000,
        40000,
        80000,
        160000,
        320000,
        640000,
        1280000,
        2560000]

    testCases.header("NumLoops", "ValueBS", "ValueMC", "TIME")

    for num_paths in num_paths_list:

        call_option = EquityDigitalOption(
            expiry_date, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type)
        value = call_option.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        start = time.time()
        value_mc = call_option.value_mc(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths)
        end = time.time()
        duration = end - start
        testCases.print(num_paths, value, value_mc, duration)

        call_option_values.append(value)
        call_option_valuesMC.append(value_mc)

#    plt.figure(figsize=(10,8))
#    plt.plot(num_paths_list, call_option_values, color = 'b', label="Call Option")
#    plt.plot(num_paths_list, call_option_valuesMC, color = 'r', label = "Call Option MC")
#    plt.xlabel("Num Loops")
#    plt.legend(loc='best')

##########################################################################

    stock_prices = range(50, 150, 50)
    call_option_values = []
    call_optionDeltas = []
    call_optionVegas = []
    call_optionThetas = []

    for stock_price in stock_prices:
        call_option = EquityDigitalOption(
            expiry_date, 100.0, OptionTypes.EUROPEAN_CALL, underlying_type)
        value = call_option.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        delta = call_option.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        vega = call_option.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        theta = call_option.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        call_option_values.append(value)
        call_optionDeltas.append(delta)
        call_optionVegas.append(vega)
        call_optionThetas.append(theta)

    put_option_values = []
    put_optionDeltas = []
    put_optionVegas = []
    put_optionThetas = []

    for stock_price in stock_prices:
        put_option = EquityDigitalOption(
            expiry_date, 100.0, OptionTypes.EUROPEAN_PUT, underlying_type)
        value = put_option.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        delta = put_option.delta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        vega = put_option.vega(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        theta = put_option.theta(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        put_option_values.append(value)
        put_optionDeltas.append(delta)
        put_optionVegas.append(vega)
        put_optionThetas.append(theta)

##########################################################################


test_EquityDigitalOption()
testCases.compareTestCases()
