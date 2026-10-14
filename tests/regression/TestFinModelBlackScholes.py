# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import matplotlib.pyplot as plt
from FinTestCases import FinTestCases, global_test_case_mode
import add_fp_to_path

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes


from financepy.utils.global_types import OptionTypes
from financepy.models.bs_american_analytic import BlackScholesAlgorithmTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.products.equity.equity_american_option import EquityAmericanOption


PLOT_GRAPHS = False


test_cases = FinTestCases(__file__, global_test_case_mode)


# TODO Complete output of results to log files

########################################################################################


def test_black_scholes():

    value_dt = Date(8, 5, 2015)
    expiry_dt = Date(15, 1, 2016)

    strike_price = 130.0
    stock_price = 127.62
    volatility = 0.20
    interest_rate = 0.001
    dividend_yield = 0.0163

    opt_type = OptionTypes.AMERICAN_CALL
    eu_option_type = OptionTypes.EUROPEAN_CALL

    # Pure American Option TREE
    amer_option = EquityAmericanOption(expiry_dt, strike_price, opt_type)

    # American with European style exercise TREE
    ameu_option = EquityAmericanOption(expiry_dt, strike_price, eu_option_type)

    # European Option and European Exercise so Black Scholes
    eu_option = EquityVanillaOption(expiry_dt, strike_price, eu_option_type)

    discount_curve = FlatDiscountCurve(
        value_dt,
        interest_rate,
        FrequencyTypes.CONTINUOUS,
    )

    dividend_curve = FlatDiscountCurve(
        value_dt,
        dividend_yield,
        FrequencyTypes.CONTINUOUS,
    )

    am_tree_value = []
    am_baw_value = []
    eu_tree_value = []
    eu_anal_value = []
    volatility = 0.20

    num_steps_per_year = range(5, 200, 1)

    test_cases.header(
        "STEPS PER YEAR",
        "AMERICAN_TREE",
        "AMERICAN_BAW",
        "EUROPEAN_TREE",
        "EUROPEAN_BS",
    )

    for num_steps in num_steps_per_year:

        model_tree = BlackScholes(volatility)
        model_anal = BlackScholes(volatility)
        model_BAW = BlackScholes(volatility)

        algo_type = BlackScholesAlgorithmTypes.CRR_TREE

        v_am = amer_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model_tree,
            algo_type=algo_type
        )

        algo_type = BlackScholesAlgorithmTypes.ANALYTICAL

        v_eu = ameu_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model_tree,
            algo_type=algo_type
        )

        algo_type = BlackScholesAlgorithmTypes.BARONE_ADESI_WHALEY

        v_bs = eu_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model_anal,
        )

        v_am_baw = amer_option.value(
            value_dt, stock_price, discount_curve, dividend_curve, model_BAW,
            algo_type=algo_type
        )

        am_tree_value.append(v_am)
        eu_tree_value.append(v_eu)
        eu_anal_value.append(v_bs)
        am_baw_value.append(v_am_baw)

        test_cases.print(num_steps, v_am, v_am_baw, v_eu, v_bs)

    if PLOT_GRAPHS:
        plt.title("American Option Price Convergence Analysis")
        plt.plot(num_steps_per_year, am_tree_value, label="American Tree")
        plt.plot(num_steps_per_year, am_baw_value, label="American BAW")
        plt.plot(num_steps_per_year, eu_tree_value, label="European Tree")
        plt.plot(num_steps_per_year, eu_anal_value, label="European Anal", lw=2)
        plt.xlabel("Num Steps")
        plt.ylabel("Value")
        plt.legend()
        plt.show()


########################################################################################

def test_barone_edesi():

    value_dt = Date(8, 5, 2015)
    expiry_dt = Date(15, 1, 2016)

    strike_price = 130.0
    stock_price = 127.62
    volatility = 0.20
    interest_rate = 0.001
    dividend_yield = 0.0163

    opt_type = OptionTypes.AMERICAN_CALL
    eu_option_type = OptionTypes.EUROPEAN_CALL

    amer_option = EquityAmericanOption(expiry_dt, strike_price, opt_type)

    ameu_option = EquityAmericanOption(expiry_dt, strike_price, eu_option_type)

    eu_option = EquityVanillaOption(expiry_dt, strike_price, eu_option_type)

    discount_curve = FlatDiscountCurve(
        value_dt, interest_rate, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )

    dividend_curve = FlatDiscountCurve(
        value_dt, dividend_yield, FrequencyTypes.CONTINUOUS, DayCountTypes.ACT_365F
    )

    model = BlackScholes(volatility)

    algo_type = BlackScholesAlgorithmTypes.DEFAULT

    v = amer_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model,
        algo_type=algo_type)

    algo_type = BlackScholesAlgorithmTypes.BARONE_ADESI_WHALEY

    v = amer_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model, algo_type=algo_type
    )

    algo_type = BlackScholesAlgorithmTypes.DEFAULT

    v = ameu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model, algo_type=algo_type
    )

    algo_type = BlackScholesAlgorithmTypes.DEFAULT

    v = eu_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model)


test_barone_edesi()
test_black_scholes()
test_cases.compare_test_cases()
