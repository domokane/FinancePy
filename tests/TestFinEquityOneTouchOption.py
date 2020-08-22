###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode

from financepy.products.equity.FinEquityOneTouchOption import FinEquityOneTouchOption
from financepy.products.equity.FinEquityOneTouchOption import FinTouchOptionPayoffTypes
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.products.equity.FinEquityModelTypes import FinEquityModelBlackScholes
from financepy.finutils.FinDate import FinDate

import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinEquityOneTouchOption():
    # Examples Haug Page 180 Table 4-22
    # Agreement not exact at t is not exactly 0.50

    valueDate = FinDate(1, 1, 2016)
    expiryDate = FinDate(2, 7, 2016)
    interestRate = 0.10
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)
    volatility = 0.20
    barrierLevel = 100.0  # H
    model = FinEquityModelBlackScholes(volatility)
    dividendYield = 0.03
    numPaths = 10000
    numStepsPerYear = 252

    stockPrice = 105.0
    paymentSize = 15.0

    print("")
    print("CASH ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_CASH_OR_NOTHING]

    for downType in downTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         downType,
                                         barrierLevel,
                                         paymentSize)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendYield,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendYield,
                              model,
                              numStepsPerYear,
                              numPaths)

        print("%60s %9.5f %9.5f" % (downType, v, v_mc))

    stockPrice = 95.0
    paymentSize = 15.0

    upTypes = [FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_HIT,
               FinTouchOptionPayoffTypes.UP_AND_IN_CASH_AT_EXPIRY,
               FinTouchOptionPayoffTypes.UP_AND_OUT_CASH_OR_NOTHING]

    for upType in upTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         upType,
                                         barrierLevel,
                                         paymentSize)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendYield,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendYield,
                              model,
                              numStepsPerYear,
                              numPaths)

        print("%60s %9.5f %9.5f" % (upType, v, v_mc))

    ###########################################################################

    stockPrice = 105.0

    print("")
    print("ASSET ONLY")

    downTypes = [FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_HIT,
                 FinTouchOptionPayoffTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
                 FinTouchOptionPayoffTypes.DOWN_AND_OUT_ASSET_OR_NOTHING]

    for downType in downTypes:

        option = FinEquityOneTouchOption(expiryDate,
                                         downType,
                                         barrierLevel)

        v = option.value(valueDate,
                         stockPrice,
                         discountCurve,
                         dividendYield,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendYield,
                              model,
                              numStepsPerYear,
                              numPaths)

        print("%60s %9.5f %9.5f" % (downType, v, v_mc))

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
                         dividendYield,
                         model)

        v_mc = option.valueMC(valueDate,
                              stockPrice,
                              discountCurve,
                              dividendYield,
                              model,
                              numStepsPerYear,
                              numPaths)

        print("%60s %9.5f %9.5f" % (upType, v, v_mc))

###############################################################################


test_FinEquityOneTouchOption()
testCases.compareTestCases()
