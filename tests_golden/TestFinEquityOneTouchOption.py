###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.products.equity.FinEquityOneTouchOption import FinEquityOneTouchOption
from financepy.products.equity.FinEquityOneTouchOption import FinTouchOptionPayoffTypes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valueDate = FinDate(1, 1, 2016)
    expiryDate = FinDate(2, 7, 2016)
    interestRate = 0.10
    volatility = 0.20
    barrierLevel = 100.0  # H
    model = FinModelBlackScholes(volatility)
    dividendYield = 0.03
    numPaths = 10000
    numStepsPerYear = 252

    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    dividendCurve = FinDiscountCurveFlat(valueDate, dividendYield)

    stockPrice = 105.0
    paymentSize = 15.0

    testCases.header("================================= CASH ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for downType in downTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         downType,
                                         barrierLevel,
                                         paymentSize)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendCurve,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stockPrice = 95.0
    paymentSize = 15.0

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_CASH_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")

    for upType in upTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         upType,
                                         barrierLevel,
                                         paymentSize)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendCurve,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    ###########################################################################

    stockPrice = 105.0

    testCases.banner("================= ASSET ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_ASSET_OR_NOTHING]

    testCases.header("TYPE", "VALUE", "VALUE_MC")
    for downType in downTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         downType,
                                         barrierLevel)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendCurve,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % downType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

    stockPrice = 95.0

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_ASSET_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_ASSET_OR_NOTHING]

    for upType in upTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         upType,
                                         barrierLevel)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendCurve,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendCurve,
                              model,
                              numStepsPerYear,
                              numPaths)

        testCases.print("%60s " % upType,
                        "%9.5f" % v,
                        "%9.5f" % v_mc)

###############################################################################


test_FinEquityOneTouchOption()
testCases.compareTestCases()
