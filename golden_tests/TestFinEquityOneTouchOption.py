###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import sys
sys.path.append("..")

from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption
from financepy.utils.global_types import TouchOptionTypes

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    value_dt = Date(1, 1, 2016)
    expiry_dt = Date(2, 7, 2016)
    interest_rate = 0.10
    volatility = 0.20
    barrier_level = 100.0  # H
    model = BlackScholes(volatility)
    dividend_yield = 0.03
    num_paths = 10000
    num_steps_per_year = 252

    discount_curve = DiscountCurveFlat(value_dt, interest_rate)
    dividend_curve = DiscountCurveFlat(value_dt, dividend_yield)

    stock_price = 105.0
    payment_size = 15.0

    test_cases.header("================================= CASH ONLY")

    downTypes = [TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT,
                 TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    test_cases.header("TYPE", "VALUE", "VALUE_MC")

    for downType in downTypes:

        option = EquityOneTouchOption(expiry_dt,
                                      downType,
                                      barrier_level,
                                      payment_size)

        v = option.value(value_dt,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(value_dt,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        test_cases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stock_price = 95.0
    payment_size = 15.0

    upTypes = [TouchOptionTypes.UP_AND_IN_CASH_AT_HIT,
               TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
               TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING]

    test_cases.header("TYPE", "VALUE", "VALUE_MC")

    for upType in upTypes:

        option = EquityOneTouchOption(expiry_dt,
                                      upType,
                                      barrier_level,
                                      payment_size)

        v = option.value(value_dt,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(value_dt,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        test_cases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    ###########################################################################

    stock_price = 105.0

    test_cases.banner("================= ASSET ONLY")

    downTypes = [TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT,
                 TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
                 TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING]

    test_cases.header("TYPE", "VALUE", "VALUE_MC")
    for downType in downTypes:

        option = EquityOneTouchOption(expiry_dt,
                                      downType,
                                      barrier_level)

        v = option.value(value_dt,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(value_dt,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        test_cases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stock_price = 95.0

    upTypes = [TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT,
               TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY,
               TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING]

    for upType in upTypes:

        option = EquityOneTouchOption(expiry_dt,
                                      upType,
                                      barrier_level)

        v = option.value(value_dt,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(value_dt,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        test_cases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

###############################################################################


test_EquityOneTouchOption()
test_cases.compareTestCases()
