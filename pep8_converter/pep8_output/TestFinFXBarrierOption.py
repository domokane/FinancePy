# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.fx.fx_barrier_option import FXBarrierOption
from financepy.products.fx.fx_barrier_option import FinFXBarrierTypes
from financepy.models.black_scholes import BlackScholes
from financepy.models.process_simulator import FinGBMNumericalScheme
from financepy.models.process_simulator import ProcessTypes


test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test__fin_fx_barrier_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    spot_fx_rate = 100.0
    currency_pair = "USDJPY"
    volatility = 0.20
    dom_interest_rate = 0.05
    for_interest_rate = 0.02
    opt_type = FinFXBarrierTypes.DOWN_AND_OUT_CALL
    notional = 100.0
    notional_currency = "USD"

    drift = dom_interest_rate - for_interest_rate
    scheme = FinGBMNumericalScheme.ANTITHETIC
    process_type = ProcessTypes.GBM
    domestic_curve = DiscountCurveFlat(value_dt, dom_interest_rate)
    foreign_curve = DiscountCurveFlat(value_dt, for_interest_rate)
    model = BlackScholes(volatility)


    import time

    start = time.time()
    num_observations_per_year = 100

    for opt_type in FinFXBarrierTypes:

        test_cases.header("Type", "k", "b", "S", "Value", "ValueMC", "TIME", "Diff")

        for spot_fx_rate in range(60, 140, 20):
            b = 110.0
            k = 100.0

            option = FXBarrierOption(
                expiry_dt,
                k,
                currency_pair,
                opt_type,
                b,
                num_observations_per_year,
                notional,
                notional_currency,
            )

            value = option.value(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            start = time.time()
            model_params = (spot_fx_rate, drift, volatility, scheme)
            value_mc = option.value_mc(
                value_dt,
                spot_fx_rate,
                dom_interest_rate,
                process_type,
                model_params,
            )

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            test_cases.print(
                opt_type,
                k,
                b,
                spot_fx_rate,
                value,
                value_mc,
                time_elapsed,
                diff,
            )

        for spot_fx_rate in range(60, 140, 20):
            b = 100.0
            k = 110.0

            option = FXBarrierOption(
                expiry_dt,
                k,
                currency_pair,
                opt_type,
                b,
                num_observations_per_year,
                notional,
                notional_currency,
            )

            value = option.value(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            start = time.time()
            model_params = (spot_fx_rate, drift, volatility, scheme)
            value_mc = option.value_mc(
                value_dt,
                spot_fx_rate,
                dom_interest_rate,
                process_type,
                model_params,
            )

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            test_cases.print(
                opt_type,
                k,
                b,
                spot_fx_rate,
                value,
                value_mc,
                time_elapsed,
                diff,
            )

    end = time.time()


    spot_fx_rates = range(50, 150, 50)
    b = 105.0

    test_cases.header("Type", "k", "b", "S:", "Value", "Delta", "Vega", "Theta")

    for opt_type in FinFXBarrierTypes:
        for spot_fx_rate in spot_fx_rates:
            barrier_option = FXBarrierOption(
                expiry_dt,
                100.0,
                currency_pair,
                opt_type,
                b,
                num_observations_per_year,
                notional,
                notional_currency,
            )

            value = barrier_option.value(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            delta = barrier_option.delta(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            vega = barrier_option.vega(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            theta = barrier_option.theta(
                value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
            )

            test_cases.print(opt_type, k, b, spot_fx_rate, value, delta, vega, theta)




################################################################################

test__fin_fx_barrier_option()
test_cases.compare_test_cases()
