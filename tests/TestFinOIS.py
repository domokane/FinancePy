###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.fin_math import ONE_MILLION
from financepy.products.rates.FinOIS import FinOIS
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.FinGlobalTypes import FinSwapTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

def test_FinFixedOIS():

    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    effective_date = Date(30, 11, 2018)
    end_date = Date(30, 11, 2023)

    end_date = effective_date.addMonths(60)
    oisRate = 0.04
    fixed_legType = FinSwapTypes.PAY
    fixedFreqType = FrequencyTypes.ANNUAL
    fixedDayCount = DayCountTypes.ACT_360
    floatFreqType = FrequencyTypes.ANNUAL
    floatDayCount = DayCountTypes.ACT_360
    floatSpread = 0.0
    notional = ONE_MILLION
    payment_lag = 1
    
    ois = FinOIS(effective_date,
                 end_date,
                 fixed_legType,
                 oisRate,
                 fixedFreqType,
                 fixedDayCount,
                 notional,
                 payment_lag,
                 floatSpread,
                 floatFreqType,
                 floatDayCount)

#    print(ois)

    valuation_date = effective_date
    marketRate = 0.05
    oisCurve = DiscountCurveFlat(valuation_date, marketRate,
                                 FrequencyTypes.ANNUAL)

    v = ois.value(effective_date, oisCurve)
    
#    print(v)
    
#    ois._fixed_leg.printValuation()
#    ois._floatLeg.printValuation()
    
    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE", v)
    
###############################################################################

test_FinFixedOIS()
testCases.compareTestCases()
