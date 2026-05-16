# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.products.fx import FXDoubleOneTouchOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.global_types import DoubleBarrierTypes

test_cases = FinTestCases(__file__, global_test_case_mode)

DEBUG_FLAG = False
########################################################################################


def test_fin_fx_double_barrier_one_touch_option():

    # Examples Haug Page 181 Table 4-23
    # Not exact due to use of dates so that t_exp ~ 0.25

    value_dt = Date(1, 1, 2016)
    expiry_dt = value_dt.add_months(3)

    spot_fx_rate = 100.0

    domestic_rate = 0.05
    foreign_rate = 0.05 - 0.03

    num_paths = 10000
    num_steps_per_year = 252 * 2

    dom_curve = DiscountCurveFlat(value_dt, domestic_rate)
    for_curve = DiscountCurveFlat(value_dt, foreign_rate)

    payment_size = 10.0

    lower_barrier_fx_rate = [80.0, 85.0, 90.0, 95.0]
    upper_barrier_fx_rate = [120.0, 115.0, 110.0, 105.0]

    for option_type in [DoubleBarrierTypes.KNOCK_OUT, DoubleBarrierTypes.KNOCK_IN]:

        for i in range(0, 4):

            k1 = lower_barrier_fx_rate[i]
            k2 = upper_barrier_fx_rate[i]

            option = FXDoubleOneTouchOption(
                expiry_dt, option_type, k1, k2, payment_size
            )

            for sigma in [0.1, 0.2, 0.3, 0.5]:

                volatility = sigma
                model = BlackScholes(volatility)

                v = option.value(value_dt, spot_fx_rate, dom_curve, for_curve, model)

                v_mc = option.value_mc(
                    value_dt,
                    spot_fx_rate,
                    dom_curve,
                    for_curve,
                    model,
                    num_steps_per_year,
                    num_paths,
                )

                #                print(option_type, k1, k2, sigma, v, v_mc)

                test_cases.header("=================================")
                test_cases.header("OPT_TYPE", "L", "U", "SIGMA", "VALUE", "VALUE_MC")
                test_cases.print(option_type, k1, k2, sigma, v, v_mc)

    #    test_cases.header("================================= CASH ONLY")
    #    test_cases.header("DELTA", "GAMMA", "VEGA")

    #        d = option.delta(value_dt, spot_fx_rate, dom_curve, for_curve, model)
    #        g = option.gamma(value_dt, spot_fx_rate, dom_curve, for_curve, model)
    #        v = option.vega(value_dt, spot_fx_rate, dom_curve, for_curve, model)

    #        test_cases.print(d, g, v)


########################################################################################

test_fin_fx_double_barrier_one_touch_option()
test_cases.compare_test_cases()
