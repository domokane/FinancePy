###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.finutils.FinMath import ONE_MILLION
from financepy.products.funding.FinFixedOIRSwap import FinFixedOIRSwap
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinGlobalTypes import FinSwapTypes

import sys
sys.path.append("..//..")


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinFixedOIRSwap():

    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    startDate = FinDate(30, 11, 2018)
    endDate = FinDate(30, 11, 2023)

    endDate = startDate.addMonths(60)
    oisRate = 0.04
    swapType = FinSwapTypes.PAYER
    fixedFreqType = FinFrequencyTypes.ANNUAL
    fixedDayCount = FinDayCountTypes.ACT_360
    floatFreqType = FinFrequencyTypes.ANNUAL
    floatDayCount = FinDayCountTypes.ACT_360
    floatSpread = 0.0
    notional = ONE_MILLION

    ois = FinFixedOIRSwap(startDate,
                          endDate,
                          swapType,
                          oisRate,
                          fixedFreqType,
                          fixedDayCount,
                          notional,
                          floatSpread,
                          floatFreqType,
                          floatDayCount)

    valueDate = FinDate(2018, 11, 30)
    marketRate = 0.05
    oisCurve = FinDiscountCurveFlat(valueDate, marketRate,
                                    FinFrequencyTypes.ANNUAL)

    v = ois.value(startDate, oisCurve)
    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE", v)
    
###############################################################################

test_FinFixedOIRSwap()
testCases.compareTestCases()
