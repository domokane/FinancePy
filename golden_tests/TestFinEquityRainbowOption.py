########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys

sys.path.append("..")

import time
from math import sqrt
import numpy as np
from financepy.products.equity.equity_rainbow_option import EquityRainbowOption
from financepy.products.equity.equity_rainbow_option import (
    EquityRainbowOptionTypes,
)
from financepy.utils.helpers import beta_vector_to_corr_matrix
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

########################################################################################


def test_EquityRainbowOption():

    #        import matplotlib.pyplot as plt

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
    num_paths_list = [10000]
    corrList = np.linspace(0.0, 0.999999, 6)
    strike = 100.0

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      CALL ON MAXIMUM")
    test_cases.banner(
        "==================================================================="
    )

    payoff_type = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

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

            v_MC = rainbow_option.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON MAX Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    ##########################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                       CALL ON MINIMUM")
    test_cases.banner(
        "==================================================================="
    )
    payoff_type = EquityRainbowOptionTypes.CALL_ON_MINIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

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

            v_MC = rainbow_option.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON MIN Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    ########################################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      PUT ON MAXIMUM")
    test_cases.banner(
        "==================================================================="
    )

    payoff_type = EquityRainbowOptionTypes.PUT_ON_MAXIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

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

            v_MC = rainbow_option.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "PUT ON MAX Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "PUT ON MAX Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    ##########################################################################

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                       PUT ON MINIMUM")
    test_cases.banner(
        "==================================================================="
    )
    payoff_type = EquityRainbowOptionTypes.PUT_ON_MINIMUM
    payoff_params = [strike]
    rainbow_option = EquityRainbowOption(
        expiry_dt, payoff_type, payoff_params, num_assets
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

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
            v_MC = rainbow_option.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "PUT ON MIN Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "PUT ON MIN Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    ##########################################################################

    num_assets = 2
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100
    strike = 100.0
    correlation = 0.50

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      CALL ON 1st")
    test_cases.banner(
        "==================================================================="
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MAXIMUM
            payoff_params1 = [strike]
            rainbowOption1 = EquityRainbowOption(
                expiry_dt, payoff_type1, payoff_params1, num_assets
            )

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [1, strike]
            rainbowOption2 = EquityRainbowOption(
                expiry_dt, payoff_type2, payoff_params2, num_assets
            )

            start = time.time()

            v = rainbowOption1.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_MC = rainbowOption2.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MAX Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON 1st Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      CALL ON 2nd")
    test_cases.banner(
        "==================================================================="
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []

    test_cases.header("NUMPATHS", "CORRELATION", "VALUE", "VALUE_MC", "TIME")

    for correlation in corrList:

        betas = np.ones(num_assets) * sqrt(correlation)
        corr_matrix = beta_vector_to_corr_matrix(betas)

        for num_paths in num_paths_list:

            payoff_type1 = EquityRainbowOptionTypes.CALL_ON_MINIMUM
            payoff_params1 = [strike]
            rainbowOption1 = EquityRainbowOption(
                expiry_dt, payoff_type1, payoff_params1, num_assets
            )

            payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
            payoff_params2 = [2, strike]
            rainbowOption2 = EquityRainbowOption(
                expiry_dt, payoff_type2, payoff_params2, num_assets
            )

            start = time.time()

            v = rainbowOption1.value(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
            )

            v_MC = rainbowOption2.value_mc(
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
            test_cases.print(num_paths, correlation, v, v_MC, duration)

            rainboxOptionValues.append(v)
            rainbowOptionValuesMC.append(v_MC)

    #    plt.figure(figsize=(10,8))
    #    plt.plot(corrList, rainboxOptionValues, color = 'r', label = "CALL ON MIN Rainbow Option Analytical")
    #    plt.plot(corrList, rainbowOptionValuesMC, 'o', color = 'b', label = "CALL ON 2nd Rainbow Option MC")
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      CALL ON 1-5")
    test_cases.banner(
        "==================================================================="
    )

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

    #    plt.figure(figsize=(10,8))

    test_cases.header(
        "NUMPATHS", "CORRELATION", "NTD", "VALUE", "VALUE_MC", "TIME"
    )

    for n in [1, 2, 3, 4, 5]:

        rainboxOptionValues = []
        rainbowOptionValuesMC = []

        payoff_type2 = EquityRainbowOptionTypes.CALL_ON_NTH
        payoff_params2 = [n, strike]
        rainbowOption2 = EquityRainbowOption(
            expiry_dt, payoff_type2, payoff_params2, num_assets
        )

        for correlation in corrList:

            betas = np.ones(num_assets) * sqrt(correlation)
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()

            v_MC = rainbowOption2.value_mc(
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
            test_cases.print(num_paths, correlation, n, v, v_MC, duration)

            rainbowOptionValuesMC.append(v_MC)

    #        plt.plot(corrList, rainbowOptionValuesMC, 'o-', label = "CALL Rainbow Option MC NTH = " + str(n))
    #    plt.xlabel("Correlation")
    #    plt.legend(loc='best')

    test_cases.banner(
        "==================================================================="
    )
    test_cases.banner("                      PUT ON 1-5")
    test_cases.banner(
        "==================================================================="
    )

    rainboxOptionValues = []
    rainbowOptionValuesMC = []
    num_paths = 10000
    num_assets = 5
    volatilities = np.ones(num_assets) * 0.3
    dividend_yields = np.ones(num_assets) * 0.01
    stock_prices = np.ones(num_assets) * 100

    #    plt.figure(figsize=(10,8))

    test_cases.header(
        "NUMPATHS", "CORRELATION", "NTD", "VALUE", "VALUE_MC", "TIME"
    )

    for n in [1, 2, 3, 4, 5]:

        rainboxOptionValues = []
        rainbowOptionValuesMC = []

        payoff_type2 = EquityRainbowOptionTypes.PUT_ON_NTH
        payoff_params2 = [n, strike]
        rainbowOption2 = EquityRainbowOption(
            expiry_dt, payoff_type2, payoff_params2, num_assets
        )

        for correlation in corrList:

            betas = np.ones(num_assets) * sqrt(correlation)
            corr_matrix = beta_vector_to_corr_matrix(betas)

            start = time.time()

            v_MC = rainbowOption2.value_mc(
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
            test_cases.print(num_paths, correlation, n, v, v_MC, duration)

            rainbowOptionValuesMC.append(v_MC)


#    plt.plot(corrList, rainbowOptionValuesMC, 'o-', label = "PUT Rainbow Option MC NTH = " + str(n))
#    plt.xlabel("Correlation")
#    plt.legend(loc='best')


########################################################################################


test_EquityRainbowOption()
test_cases.compareTestCases()
