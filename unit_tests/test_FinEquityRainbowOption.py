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

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
interest_rate = 0.05

discount_curve = DiscountCurveFlat(value_dt, interest_rate)

num_assets = 2
volatilities = np.ones(num_assets) * 0.3

dividend_yields = np.ones(num_assets) * 0.01

dividend_curves = []
for q in dividend_yields:
    dividend_curve = DiscountCurveFlat(value_dt, q)
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
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    v_MC = rainbowOption.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    assert round(v, 4) == 21.4034
    assert round(v_MC, 4) == 21.3580


def test_call_on_min():
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    v_MC = rainbowOption.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    assert round(v, 4) == 5.7817
    assert round(v_MC, 4) == 5.8022


def test_put_on_max():
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    v_MC = rainbowOption.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    assert round(v, 4) == 4.6493
    assert round(v_MC, 4) == 4.6484


def test_put_on_min():
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoff_params = [strike]
    rainbowOption = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    betas = np.ones(num_assets) * sqrt(correlation)
    corr_matrix = beta_vector_to_corr_matrix(betas)

    v = rainbowOption.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    v_MC = rainbowOption.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    assert round(v, 4) == 14.8750
    assert round(v_MC, 4) == 14.7673


def test_call_on_nth():
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    payoff_type = EquityRainbowOptionTypes.CALL_ON_NTH

    expected_results = [34.4109, 18.1054, 9.7345, 4.7557, 1.6387]

    for n in [1, 2, 3, 4, 5]:
        print(n)
        payoff_params = [n, strike]
        rainbowOption = EquityRainbowOption(
            expiry_dt, payoff_type, payoff_params, num_assets
        )

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        v_MC = rainbowOption.value_mc(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths,
        )

        assert round(v_MC, 4) == expected_results[n - 1]


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
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    payoff_type = EquityRainbowOptionTypes.PUT_ON_NTH

    expected_results = [1.4413, 3.9902, 7.7395, 13.3805, 22.4640]

    for n in [1, 2, 3, 4, 5]:
        print(n)
        payoff_params = [n, strike]
        rainbowOption = EquityRainbowOption(
            expiry_dt, payoff_type, payoff_params, num_assets
        )

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        v_MC = rainbowOption.value_mc(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths,
        )

        assert round(v_MC, 4) == expected_results[n - 1]
