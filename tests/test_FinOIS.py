###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.products.rates.ois import OIS
from financepy.utils.math import ONE_MILLION


def test_FinFixedOIS():
    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    effective_date = Date(30, 11, 2018)
    end_date = Date(30, 11, 2023)

    end_date = effective_date.add_months(60)
    oisRate = 0.04
    fixed_leg_type = SwapTypes.PAY
    fixedFreqType = FrequencyTypes.ANNUAL
    fixedDayCount = DayCountTypes.ACT_360
    floatFreqType = FrequencyTypes.ANNUAL
    floatDayCount = DayCountTypes.ACT_360
    float_spread = 0.0
    notional = ONE_MILLION
    payment_lag = 1

    ois = OIS(effective_date,
              end_date,
              fixed_leg_type,
              oisRate,
              fixedFreqType,
              fixedDayCount,
              notional,
              payment_lag,
              float_spread,
              floatFreqType,
              floatDayCount)

    valuation_date = effective_date
    marketRate = 0.05
    oisCurve = DiscountCurveFlat(valuation_date, marketRate,
                                 FrequencyTypes.ANNUAL)

    v = ois.value(effective_date, oisCurve)
    assert round(v, 4) == 43915.6019
