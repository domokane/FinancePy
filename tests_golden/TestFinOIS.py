###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.utils.math import ONE_MILLION
from financepy.products.rates.ois import OIS
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from FinTestCases import FinTestCases, globalTestCaseMode


testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFixedOIS():

    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    effective_date = Date(30, 11, 2018)
    end_date = Date(30, 11, 2023)

    end_date = effective_date.add_months(60)
    oisRate = 0.04
    fixed_leg_type = SwapTypes.PAY
    fixed_freq_type = FrequencyTypes.ANNUAL
    fixedDayCount = DayCountTypes.ACT_360
    float_freq_type = FrequencyTypes.ANNUAL
    floatDayCount = DayCountTypes.ACT_360
    float_spread = 0.0
    notional = ONE_MILLION
    payment_lag = 1

    ois = OIS(effective_date,
              end_date,
              fixed_leg_type,
              oisRate,
              fixed_freq_type,
              fixedDayCount,
              notional,
              payment_lag,
              float_spread,
              float_freq_type,
              floatDayCount)

#    print(ois)

    value_date = effective_date
    marketRate = 0.05
    oisCurve = DiscountCurveFlat(value_date,
                                 marketRate,
                                 FrequencyTypes.ANNUAL)

    v = ois.value(effective_date, oisCurve)

#    print(v)

#    ois._fixed_leg.print_valuation()
#    ois._float_leg.print_valuation()

    testCases.header("LABEL", "VALUE")
    testCases.print("SWAP_VALUE", v)

###############################################################################


test_FinFixedOIS()
testCases.compareTestCases()
