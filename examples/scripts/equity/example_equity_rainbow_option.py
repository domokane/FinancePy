# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path

_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause

_install_double_click_pause()
import numpy as np
import time
from math import sqrt

import add_fp_to_path

from financepy.products.equity.equity_rainbow_option import EquityRainbowOption
from financepy.products.equity.equity_rainbow_option import (
    EquityRainbowOptionTypes,
)
from financepy.utils.helpers import beta_vector_to_corr_matrix
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date

########################################################################################


def test_equity_rainbow_option():

    #        import matplotlib.pyplot as plt

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    interest_rate = 0.05

    discount_curve = FlatDiscountCurve(value_dt, interest_rate)

    num_assets = 2
    volatilities = np.ones(num_assets) * 0.3

    dividend_yields = np.ones(num_assets) * 0.01

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = FlatDiscountCurve(value_dt, q)
        dividend_curves.append(dividend_curve)

    stock_prices = np.ones(num_assets) * 100
    num_paths_list = [10000]
    corr_list = np.linspace(0.0, 0.999999, 6)
    strike = 100.0

    print("===================================================================")
    print("                      CALL ON MAXIMUM")
    print("===================================================================")

    payoff_type = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(expiry_dt, payoff_type, payoff_params, num_assets)

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            start = time.time()
            v = rainbow_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = rainbow_option.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "CALL ON MAX Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                       CALL ON MINIMUM")
    print("===================================================================")
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(expiry_dt, payoff_type, payoff_params, num_assets)

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            start = time.time()

            v = rainbow_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = rainbow_option.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "CALL ON MIN Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                      PUT ON MAXIMUM")
    print("===================================================================")

    payoff_type = EquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(expiry_dt, payoff_type, payoff_params, num_assets)

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            start = time.time()

            v = rainbow_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = rainbow_option.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "PUT ON MAX Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "PUT ON MAX Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                       PUT ON MINIMUM")
    print("===================================================================")
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(expiry_dt, payoff_type, payoff_params, num_assets)

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            start = time.time()
            v = rainbow_option.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )
            v_mc = rainbow_option.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "PUT ON MIN Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "PUT ON MIN Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    num_assets = 2
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100
    strike = 100.0
    correlation = 0.50

    print("===================================================================")
    print("                      CALL ON 1st")
    print("===================================================================")

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
            payoff_params1 = [strike]
            rainbow_option1 = EquityRainbowOption(expiry_dt, payoff_type1, payoff_params1, num_assets)

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [1, strike]
            rainbow_option2 = EquityRainbowOption(expiry_dt, payoff_type2, payoff_params2, num_assets)

            start = time.time()

            v = rainbow_option1.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = rainbow_option2.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "CALL ON 1st Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                      CALL ON 2nd")
    print("===================================================================")

    rainbox_option_values = []
    rainbow_option_values_mc = []

    print("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corr_list:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MINIMUM
            payoff_params1 = [strike]
            rainbow_option1 = EquityRainbowOption(expiry_dt, payoff_type1, payoff_params1, num_assets)

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [2, strike]
            rainbow_option2 = EquityRainbowOption(expiry_dt, payoff_type2, payoff_params2, num_assets)

            start = time.time()

            v = rainbow_option1.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_mc = rainbow_option2.value_mc(
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
            print(num_paths, correlation, v, v_mc, duration)

            rainbox_option_values.append(v)
            rainbow_option_values_mc.append(v_mc)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corr_list, rainbox_option_values, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
    #    plt.plot(corr_list, rainbow_option_values_mc, 'o', color = 'b', label = "CALL ON 2nd Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                      CALL ON 1-5")
    print("===================================================================")

    rainbox_option_values = []
    rainbow_option_values_mc = []
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    dividend_curves = []
    for q in dividend_yields:
        dividend_curve = FlatDiscountCurve(value_dt, q)
        dividend_curves.append(dividend_curve)

    #    plt.figure(figsize=(10,8))

    print("NUMPATHS", "CORRELATION", "NTD", "VALUE", "VALUE_MC", "TIME")

    for n in [1, 2, 3, 4, 5]:

        rainbox_option_values = []
        rainbow_option_values_mc = []

        payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
        payoff_params2 = [n, strike]
        rainbow_option2 = EquityRainbowOption(expiry_dt, payoff_type2, payoff_params2, num_assets)

        for correlation in corr_list:

            betas = np.ones(num_assets) * sqrt(correlation)
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()

            v_mc = rainbow_option2.value_mc(
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
            print(num_paths, correlation, n, v, v_mc, duration)

            rainbow_option_values_mc.append(v_mc)

    #        plt.plot(corr_list, rainbow_option_values_mc, 'o-', label = "CALL Rainbow Option MC NTH = " + str(n))
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    print("===================================================================")
    print("                      PUT ON 1-5")
    print("===================================================================")

    rainbox_option_values = []
    rainbow_option_values_mc = []
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    #    plt.figure(figsize=(10,8))

    print("NUMPATHS", "CORRELATION", "NTD", "VALUE", "VALUE_MC", "TIME")

    for n in [1, 2, 3, 4, 5]:

        rainbox_option_values = []
        rainbow_option_values_mc = []

        payoff_type2 = EquityRainbowOptionTypes.PUT_ON_NTH
        payoff_params2 = [n, strike]
        rainbow_option2 = EquityRainbowOption(expiry_dt, payoff_type2, payoff_params2, num_assets)

        for correlation in corr_list:

            betas = np.ones(num_assets) * sqrt(correlation)
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()

            v_mc = rainbow_option2.value_mc(
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
            print(num_paths, correlation, n, v, v_mc, duration)

            rainbow_option_values_mc.append(v_mc)


########################################################################################

#    plt.plot(corr_list, rainbow_option_values_mc, 'o-', label = "PUT Rainbow Option MC NTH = " + str(n))
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')


test_equity_rainbow_option()
