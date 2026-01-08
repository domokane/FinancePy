# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
import time

import add_fp_to_path

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.models.process_simulator import FinGBMNumericalScheme

from financepy.products.equity.equity_barrier_option import BarrierTypes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_equity_barrier_option():

    value_dt = Date(1, 1, 2015)
    expiry_dt = Date(1, 1, 2016)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    opt_type = BarrierTypes.DOWN_AND_OUT_CALL

    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)
    model = BlackScholes(volatility)

    start = time.time()
    num_obs_per_year = 252
    num_paths = 1000
    seed = 42

    test_cases.header("Type", "K", "B", "S:", "Value:", "ValueMC", "Diff", "TIME")

    for opt_type in BarrierTypes:
        for stock_price in [80, 100, 120]:

            barrier = 110.0
            strike = 100.0

            option = EquityBarrierOption(
                expiry_dt, strike, opt_type, barrier, num_obs_per_year
            )

            value = option.value(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )

            start = time.time()

            value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                num_obs_per_year=num_obs_per_year,
                num_paths=num_paths,
                seed=seed,
            )

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = value_mc - value

            test_cases.print(
                opt_type,
                strike,
                barrier,
                stock_price,
                value,
                value_mc,
                diff,
                time_elapsed,
            )

        for stock_price in [80, 100, 120]:

            strike = 110.0
            barrier = 100.0

            option = EquityBarrierOption(
                expiry_dt, strike, opt_type, barrier, num_obs_per_year
            )

            value = option.value(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )
            start = time.time()

            test_value_mc = option.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                num_obs_per_year=num_obs_per_year,
                num_paths=num_paths,
                seed=seed,
            )

            end = time.time()
            time_elapsed = round(end - start, 3)
            diff = test_value_mc - value

            test_cases.print(
                opt_type,
                strike,
                barrier,
                stock_price,
                value,
                test_value_mc,
                diff,
                time_elapsed,
            )

        end = time.time()

    stock_prices = [80, 100, 120]
    barrier = 105.0
    strike = 100.0

    test_cases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    for opt_type in BarrierTypes:

        for stock_price in stock_prices:

            barrier_option = EquityBarrierOption(
                expiry_dt, strike, opt_type, barrier, num_obs_per_year
            )

            value = barrier_option.value(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )
            delta = barrier_option.delta(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )
            vega = barrier_option.vega(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )
            theta = barrier_option.theta(
                value_dt, stock_price, discount_curve, dividend_curve, model
            )

            test_cases.print(
                opt_type, strike, barrier, stock_price, value, delta, vega, theta
            )

    stock_prices = [80, 100, 120]
    barrier = 105.0

    test_cases.header("Type", "K", "B", "S:", "Value", "Delta", "Vega", "Theta")

    barrier_option = EquityBarrierOption(
        expiry_dt, 100.0, opt_type, barrier, num_obs_per_year
    )


########################################################################################

test_equity_barrier_option()
test_cases.compare_test_cases()
