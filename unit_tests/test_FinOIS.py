########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

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

    ois = OIS(
        effective_dt,
        end_dt,
        fixed_leg_type,
        oisRate,
        fixed_freq_type,
        fixedDayCount,
        notional,
        payment_lag,
        float_spread,
        float_freq_type,
        floatDayCount,
    )

    value_dt = effective_dt
    marketRate = 0.05
    ois_curve = DiscountCurveFlat(value_dt, marketRate, FrequencyTypes.ANNUAL)

    v = ois.value(effective_dt, ois_curve)
    assert round(v, 4) == 43915.6019
