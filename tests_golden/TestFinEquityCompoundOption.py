###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_compound_option import EquityCompoundOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_EquityCompoundOption():

    valuation_date = Date(1, 1, 2015)
    expiry_date1 = Date(1, 1, 2017)
    expiry_date2 = Date(1, 1, 2018)
    k1 = 5.0
    k2 = 95.0
    stock_price = 85.0
    volatility = 0.15
    interest_rate = 0.035
    dividend_yield = 0.01

    model = BlackScholes(volatility)
    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    num_steps_list = [100, 200, 500, 1000, 2000, 5000]

    ###########################################################################

    stock_price = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S",
                     "TreeSteps", "Exact", "TreeValue")

    for option_type1 in [
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT]:
        for option_type2 in [
                OptionTypes.EUROPEAN_CALL,
                OptionTypes.EUROPEAN_PUT]:

            cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                              expiry_date2, option_type2, k2)

            for num_steps in num_steps_list:

                value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                                         dividend_curve, model)

                values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                                dividend_curve, model, num_steps)

                testCases.print(option_type1, option_type2, k1, k2, stock_price,
                                num_steps, value, values[0])

    ###########################################################################

    stock_price = 85.0

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S",
                     "TreeSteps", "Exact", "TreeValue")

    for option_type1 in [
            OptionTypes.AMERICAN_CALL,
            OptionTypes.AMERICAN_PUT]:
        for option_type2 in [
                OptionTypes.AMERICAN_CALL,
                OptionTypes.AMERICAN_PUT]:

            cmpdOption = EquityCompoundOption(expiry_date1, option_type1, k1,
                                              expiry_date2, option_type2, k2)

            for num_steps in num_steps_list:

                value = cmpdOption.value(valuation_date, stock_price, discount_curve,
                                         dividend_curve, model, num_steps)

                values = cmpdOption._value_tree(valuation_date, stock_price, discount_curve,
                                                dividend_curve, model, num_steps)

                testCases.print(option_type1, option_type2, k1, k2, stock_price,
                                num_steps, value, values[0])

    ###########################################################################

    testCases.header("TYPE1", "TYPE2", "K1", "K2", "S", "Exact", "TreeSteps",
                     "TreeValue", "Diff", "DELTA", "GAMMA", "THETA")

    for option_type1 in [
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT]:
        for option_type2 in [
                OptionTypes.EUROPEAN_CALL,
                OptionTypes.EUROPEAN_PUT]:

            cmpdOption = EquityCompoundOption(
                expiry_date1, option_type1, k1,
                expiry_date2, option_type2, k2)
            stock_prices = range(70, 100, 10)

            for stock_price in stock_prices:
                value = cmpdOption.value(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model)
                delta = cmpdOption.delta(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model)
                vega = cmpdOption.vega(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model)
                theta = cmpdOption.theta(
                    valuation_date,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model)

                values = cmpdOption._value_tree(valuation_date, stock_price,
                                                discount_curve, dividend_curve,
                                                model)

                diff = value - values[0]

                testCases.print(
                    option_type1,
                    option_type2,
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


test_EquityCompoundOption()
testCases.compareTestCases()
