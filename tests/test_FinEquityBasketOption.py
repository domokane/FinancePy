###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.utils.helpers import beta_vector_to_corr_matrix
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.global_types import OptionTypes
from financepy.products.equity.equity_basket_option import EquityBasketOption
import numpy as np


valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
volatility = 0.30
interest_rate = 0.05
discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
num_assets = 5
beta = 0.999999
betas = np.ones(num_assets) * beta
corr_matrix = beta_vector_to_corr_matrix(betas)
num_paths = 10000


def test_homogeneous_call():
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    call_option = EquityBasketOption(
        expiry_date, 100.0, OptionTypes.EUROPEAN_CALL, num_assets)
    value = call_option.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    assert round(value, 4) == 13.6164

    value_mc = call_option.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(value_mc, 4) == 13.5338


def test_homogeneous_put():
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    call_option = EquityBasketOption(
        expiry_date, 100.0, OptionTypes.EUROPEAN_PUT, num_assets)
    value = call_option.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    assert round(value, 4) == 9.7344

    value_mc = call_option.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(value_mc, 4) == 9.6986


def test_inhomogeneous_call():
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    call_option = EquityBasketOption(
        expiry_date, 100.0, OptionTypes.EUROPEAN_CALL, num_assets)
    value = call_option.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    assert round(value, 4) == 13.6783

    value_mc = call_option.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(value_mc, 4) == 13.5460


def test_inhomogeneous_put():
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    call_option = EquityBasketOption(
        expiry_date, 100.0, OptionTypes.EUROPEAN_PUT, num_assets)
    value = call_option.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    assert round(value, 4) == 7.9126

    value_mc = call_option.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(value_mc, 4) == 7.8216
