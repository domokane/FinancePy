###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption
from financepy.products.equity.equity_one_touch_option import FinTouchOptionPayoffTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_EquityOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(2, 7, 2016)
    interest_rate = 0.10
    volatility = 0.20
    barrier_level = 100.0  # H
    model = BlackScholes(volatility)
    dividend_yield = 0.03
    num_paths = 10000
    num_steps_per_year = 252

    discount_curve = DiscountCurveFlat(valuation_date, interest_rate)
    dividend_curve = DiscountCurveFlat(valuation_date, dividend_yield)

    stock_price = 105.0
    payment_size = 15.0

    testCases.header("================================= CASH ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for downType in downTypes:

        option = EquityOneTouchOption(expiry_date,
                                      downType,
                                      barrier_level,
                                      payment_size)

        v = option.value(valuation_date,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stock_price = 95.0
    payment_size = 15.0

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for upType in upTypes:

        option = EquityOneTouchOption(expiry_date,
                                      upType,
                                      barrier_level,
                                      payment_size)

        v = option.value(valuation_date,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    ###########################################################################

    stock_price = 105.0

    testCases.banner("================= ASSET ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_ASSET_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")
    for downType in downTypes:

        option = EquityOneTouchOption(expiry_date,
                                      downType,
                                      barrier_level)

        v = option.value(valuation_date,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stock_price = 95.0

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_ASSET_OR_NOTHING]

    for upType in upTypes:

        option = EquityOneTouchOption(expiry_date,
                                      upType,
                                      barrier_level)

        v = option.value(valuation_date,
                         stock_price,
                         discount_curve,
                         dividend_curve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               stock_price,
                               discount_curve,
                               dividend_curve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

###############################################################################


test_EquityOneTouchOption()
testCases.compareTestCases()
