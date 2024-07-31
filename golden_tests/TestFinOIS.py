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


test_cases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinFixedOIS():

    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    effective_dt = Date(30, 11, 2018)
    end_dt = Date(30, 11, 2023)

    end_dt = effective_dt.add_months(60)
    oisRate = 0.04
    fixed_leg_type = SwapTypes.PAY
    fixed_freq_type = FrequencyTypes.ANNUAL
    fixedDayCount = DayCountTypes.ACT_360
    float_freq_type = FrequencyTypes.ANNUAL
    floatDayCount = DayCountTypes.ACT_360
    float_spread = 0.0
    notional = ONE_MILLION
    payment_lag = 1

    ois = OIS(effective_dt,
              end_dt,
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

    value_dt = effective_dt
    marketRate = 0.05
    oisCurve = DiscountCurveFlat(value_dt,
                                 marketRate,
                                 FrequencyTypes.ANNUAL)

    v = ois.value(effective_dt, oisCurve)

#    print(v)

#    ois._fixed_leg.print_valuation()
#    ois._float_leg.print_valuation()

    test_cases.header("LABEL", "VALUE")
    test_cases.print("SWAP_VALUE", v)

###############################################################################


test_FinFixedOIS()
test_cases.compareTestCases()
