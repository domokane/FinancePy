###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from financepy.products.equity.equity_index_option import EquityIndexOption
from financepy.utils.global_types import OptionTypes
from financepy.models.black import Black, BlackTypes
from financepy.models.equity_crr_tree import crr_tree_val_avg
from financepy.utils.global_vars import gDaysInYear


def test_equity_european_index_option_price():
    # https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1540-6261.1986.tb04495.x
    # exiray_date is set so that time to maturity becomes 0.25
    valuation_date = Date(8, 5, 2015)
    expiry_date = Date(7, 8, 2015, hh=6)
    discountRate = 0.08
    volatility = 0.15
    discount_curve = DiscountCurveFlat(valuation_date, discountRate)
    strike = 100.0
    future_prices = [80.0, 90.0, 100.0, 110.0, 120.0]
    # expected result
    expected_european_call_prices = [
        0.0027, 0.2529, 2.9321, 10.1752, 19.6239]
    expected_european_put_prices = [
        19.6067, 10.0549, 2.9321, 0.3732, 0.0199]

    # costruct each opton
    EuropeanCallOption = EquityIndexOption(
        expiry_date, strike, OptionTypes.EUROPEAN_CALL)
    EuropeanPutOption = EquityIndexOption(
        expiry_date, strike, OptionTypes.EUROPEAN_PUT)

    # construct black model
    modelBsAnalytical = Black(
        volatility, implementation_type=BlackTypes.ANALYTICAL)
    for i in range(len(future_prices)):
        price = EuropeanCallOption.value(
            valuation_date, future_prices[i], discount_curve, modelBsAnalytical)
        assert round(price, 1) == round(expected_european_call_prices[i], 1)
        price = EuropeanPutOption.value(
            valuation_date, future_prices[i], discount_curve, modelBsAnalytical)
        assert round(price, 1) == round(expected_european_put_prices[i], 1)


def test_equity_american_index_option_price():
    # Just check crr_tree_val_avg can called correctly from EquityIndexOption
    # in case of American exercise
    valuation_date = Date(8, 5, 2015)
    expiry_date = Date(7, 8, 2015, hh=6)
    time_to_expiry = (expiry_date-valuation_date) / gDaysInYear
    discountRate = 0.08
    volatility = 0.15
    discount_curve = DiscountCurveFlat(valuation_date, discountRate)
    strike = 100.0
    future_prices = [80.0, 90.0, 100.0, 110.0, 120.0]
    num_steps = 200
    # costruct each opton
    AmericanCallOption = EquityIndexOption(
        expiry_date, strike, OptionTypes.AMERICAN_CALL)
    AmericanPutOption = EquityIndexOption(
        expiry_date, strike, OptionTypes.AMERICAN_PUT)
    # construct black model
    modelBsCrrTree = Black(
        volatility, implementation_type=BlackTypes.CRR_TREE, num_steps=num_steps)
    for i in range(len(future_prices)):
        price = AmericanCallOption.value(
            valuation_date, future_prices[i], discount_curve, modelBsCrrTree)
        result = crr_tree_val_avg(
            future_prices[i], 0.0, 0.0, volatility, num_steps,
            time_to_expiry, OptionTypes.AMERICAN_CALL.value, strike)
        price_expected = result['value']
        assert round(price, 4) == round(price_expected, 4)
        price = AmericanPutOption.value(
            valuation_date, future_prices[i], discount_curve, modelBsCrrTree)
        result = crr_tree_val_avg(
            future_prices[i], 0.0, 0.0, volatility, num_steps,
            time_to_expiry, OptionTypes.AMERICAN_PUT.value, strike)
        price_expected = result['value']
        assert round(price, 4) == round(price_expected, 4)
