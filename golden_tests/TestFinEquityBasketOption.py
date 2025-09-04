# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.products.equity.equity_basket_option import EquityBasketOption
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.helpers import beta_vector_to_corr_matrix
from financepy.utils.date import Date

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_basket_option():

    import time

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    volatility = 0.30
    interest_rate = 0.05
    discount_curve = DiscountCurveFlat(value_dt, interest_rate)

    # Homogeneous Basket

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    beta_list = np.linspace(0.0, 0.999999, 11)

    test_cases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in beta_list:
        for num_paths in [10000]:
            call_option = EquityBasketOption(
                expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, num_assets
            )
            betas = np.ones(num_assets) * beta
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()
            v = call_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = call_option.value_mc(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
                num_paths,
            )
            end = time.time()
            duration = end - start
            test_cases.print(num_paths, beta, v, v_mc, duration)

    # INHomogeneous Basket

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    beta_list = np.linspace(0.0, 0.999999, 11)

    test_cases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in beta_list:

        for num_paths in [10000]:

            call_option = EquityBasketOption(
                expiry_dt, 100.0, OptionTypes.EUROPEAN_CALL, num_assets
            )
            betas = np.ones(num_assets) * beta
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()

            v = call_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = call_option.value_mc(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
                num_paths,
            )

            end = time.time()
            duration = end - start
            test_cases.print(num_paths, beta, v, v_mc, duration)

    # Homogeneous Basket

    num_assets = 5
    volatilities = np.ones(num_assets) * volatility
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100
    beta_list = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    test_cases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in beta_list:
        for num_paths in [10000]:
            call_option = EquityBasketOption(
                expiry_dt, 100.0, OptionTypes.EUROPEAN_PUT, num_assets
            )
            betas = np.ones(num_assets) * beta
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()
            v = call_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )
            v_mc = call_option.value_mc(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
                num_paths,
            )
            end = time.time()
            duration = end - start
            test_cases.print(num_paths, beta, v, v_mc, duration)

    # INHomogeneous Basket

    num_assets = 5
    volatilities = np.array([0.3, 0.2, 0.25, 0.22, 0.4])
    dividend_yields = np.array([0.01, 0.02, 0.04, 0.01, 0.02])
    stock_prices = np.array([100, 105, 120, 100, 90])
    beta_list = np.linspace(0.0, 0.999999, 11)

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = DiscountCurveFlat(value_dt, q)
        dividend_curves.append(dividend_curve)

    test_cases.header("NumPaths", "Beta", "Value", "ValueMC", "TIME")

    for beta in beta_list:

        for num_paths in [10000]:

            call_option = EquityBasketOption(
                expiry_dt, 100.0, OptionTypes.EUROPEAN_PUT, num_assets
            )
            betas = np.ones(num_assets) * beta
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()
            v = call_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )
            v_mc = call_option.value_mc(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
                num_paths,
            )
            end = time.time()
            duration = end - start
            test_cases.print(num_paths, beta, v, v_mc, duration)


########################################################################################

test_equity_basket_option()
test_cases.compare_test_cases()
