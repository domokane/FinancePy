###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption
from financepy.products.equity.equity_one_touch_option import FinTouchOptionPayoffTypes
from financepy.market.discount.curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFXOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valueDate = FinDate(1, 1, 2016)
    expiryDate = FinDate(2, 7, 2016)
    volatility = 0.20
    barrier_level = 1.0  # H
    model = BlackScholes(volatility)

    domesticRate = 0.10
    foreignRate = 0.03

    numPaths = 20000
    numStepsPerYear = 252

    domCurve = FinDiscountCurveFlat(valueDate, domesticRate)
    forCurve = FinDiscountCurveFlat(valueDate, foreignRate)

    spot_fx_rate = 1.050
    paymentSize = 1.5

    testCases.header("================================= CASH ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for downType in downTypes:

        option = EquityOneTouchOption(expiryDate,
                                         downType,
                                         barrier_level,
                                         paymentSize)

        v = option.value(valueDate,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valueDate,
                              spot_fx_rate,
                              domCurve,
                              forCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    spot_fx_rate = 0.950
    paymentSize = 1.5

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for upType in upTypes:

        option = EquityOneTouchOption(expiryDate,
                                         upType,
                                         barrier_level,
                                         paymentSize)

        v = option.value(valueDate,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valueDate,
                              spot_fx_rate,
                              domCurve,
                              forCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

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

        option = EquityOneTouchOption(expiryDate,
                                         downType,
                                         barrier_level)

        v = option.value(valueDate,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valueDate,
                              spot_fx_rate,
                              domCurve,
                              forCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    spot_fx_rate = 0.950

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_ASSET_OR_NOTHING]

    for upType in upTypes:

        option = EquityOneTouchOption(expiryDate,
                                         upType,
                                         barrier_level)

        v = option.value(valueDate,
                         spot_fx_rate,
                         domCurve,
                         forCurve,
                         model)

        v_mc = option.value_mc(valueDate,
                              spot_fx_rate,
                              domCurve,
                              forCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

###############################################################################


test_FinFXOneTouchOption()
testCases.compareTestCases()
