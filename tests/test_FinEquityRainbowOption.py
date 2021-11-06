###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.helpers import beta_vector_to_corr_matrix
from financepy.products.equity.equity_rainbow_option import EquityRainbowOptionTypes
from financepy.products.equity.equity_rainbow_option import EquityRainbowOption
import numpy as np
from math import sqrt

valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 1, 2016)
interest_rate = 0.05

discount_curve = DiscountCurveFlat(valuation_date, interest_rate)

num_assets = 2
volatilities = np.ones(num_assets) * 0.3

dividend_yields = np.ones(num_assets) * 0.01

dividend_curves = []
for q in dividend_yields:
    dividend_curve = DiscountCurveFlat(valuation_date, q)
    dividend_curves.append(dividend_curve)

stock_prices = np.ones(num_assets) * 100
num_paths = 10000
corrList = np.linspace(0.0, 0.999999, 6)
strike = 100.0

correlation = 0.39999960


def test_call_on_max():
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_date, payoff_type, payoff_params, num_assets)

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    v_MC = rainbowOption.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(v, 4) == 21.4034
    assert round(v_MC, 4) == 21.5586


def test_call_on_min():
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_date, payoff_type, payoff_params, num_assets)

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    v_MC = rainbowOption.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(v, 4) == 5.7817
    assert round(v_MC, 4) == 5.8795


def test_put_on_max():
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_date, payoff_type, payoff_params, num_assets)

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    v_MC = rainbowOption.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(v, 4) == 4.6493
    assert round(v_MC, 4) == 4.6839


def test_put_on_min():
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_date, payoff_type, payoff_params, num_assets)

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix)

    v_MC = rainbowOption.value_mc(
        valuation_date,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths)

    assert round(v, 4) == 14.8750
    assert round(v_MC, 4) == 14.8747


def test_call_on_nth():
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    payoff_type = EquityRainbowOptionTypes.CALL_ON_NTH

    expected_results = [
        34.1760,
        17.8990,
        9.6716,
        4.7205,
        1.6182
    ]

    for n in [1, 2, 3, 4, 5]:
        print(n)
        payoff_params = [n, strike]
        rainbowOption = EquityRainbowOption(
            expiry_date, payoff_type, payoff_params, num_assets)

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        v_MC = rainbowOption.value_mc(
            valuation_date,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths)

        assert round(v_MC, 4) == expected_results[n-1]


def test_put_on_nth():
    rainboxOptionValues = []
    rainbowOptionValuesMC = []
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(valuation_date, q)
        dividend_curves.append(dividend_curve)

    payoff_type = EquityRainbowOptionTypes.PUT_ON_NTH

    expected_results = [
        1.4277,
        3.9644,
        7.6921,
        13.2466,
        22.3467
    ]

    for n in [1, 2, 3, 4, 5]:
        print(n)
        payoff_params = [n, strike]
        rainbowOption = EquityRainbowOption(
            expiry_date, payoff_type, payoff_params, num_assets)

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        v_MC = rainbowOption.value_mc(
            valuation_date,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths)

        assert round(v_MC, 4) == expected_results[n-1]
