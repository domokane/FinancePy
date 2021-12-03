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


def test_FinFXOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valuation_date = Date(1, 1, 2016)
    expiry_date = Date(2, 7, 2016)
    volatility = 0.20
    barrier_level = 1.0  # H
    model = BlackScholes(volatility)

    domesticRate = 0.10
    foreignRate = 0.03

    num_paths = 20000
    num_steps_per_year = 252

    domCurve = DiscountCurveFlat(valuation_date, domesticRate)
    forCurve = DiscountCurveFlat(valuation_date, foreignRate)

    spot_fx_rate = 1.050
    payment_size = 1.5

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
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    spot_fx_rate = 0.950
    payment_size = 1.5

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
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    ###########################################################################

    spot_fx_rate = 1.050

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
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    spot_fx_rate = 0.950

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_ASSET_OR_NOTHING]

    for upType in upTypes:

        option = EquityOneTouchOption(expiry_date,
                                      upType,
                                      barrier_level)

        v = option.value(valuation_date,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valuation_date,
                               spot_fx_rate,
                               domCurve,
                               forCurve,
                               model,
                               num_steps_per_year,
                               num_paths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

###############################################################################


test_FinFXOneTouchOption()
testCases.compareTestCases()
