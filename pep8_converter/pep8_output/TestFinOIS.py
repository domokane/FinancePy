# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

from financepy.utils.math import ONE_MILLION
from financepy.products.rates.ois import OIS
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
from FinTestCases import FinTestCases, global_test_case_mode


test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test__fin_fixed_ois():

    # Here I follow the example in
    # https://blog.deriscope.com/index.php/en/excel-quantlib-overnight-index-swap

    effective_dt = Date(30, 11, 2018)
    end_dt = Date(30, 11, 2023)

    end_dt = effective_dt.add_months(60)
    ois_rate = 0.04
    fixed_leg_type = SwapTypes.PAY
    fixed_freq_type = FrequencyTypes.ANNUAL
    fixed_day_count = DayCountTypes.ACT_360
    float_freq_type = FrequencyTypes.ANNUAL
    float_day_count = DayCountTypes.ACT_360
    float_spread = 0.0
    notional = ONE_MILLION
    payment_lag = 1

    ois = OIS(
        effective_dt,
        end_dt,
        fixed_leg_type,
        ois_rate,
        fixed_freq_type,
        fixed_day_count,
        notional,
        payment_lag,
        float_spread,
        float_freq_type,
        float_day_count,
    )

    #    print(ois)

    value_dt = effective_dt
    market_rate = 0.05
    ois_curve = DiscountCurveFlat(value_dt, market_rate, FrequencyTypes.ANNUAL)

    v = ois.value(effective_dt, ois_curve)

    #    print(v)

    #    ois._fixed_leg.print_valuation()
    #    ois._float_leg.print_valuation()

    test_cases.header("LABEL", "VALUE")
    test_cases.print("SWAP_VALUE", v)




################################################################################

test__fin_fixed_ois()
test_cases.compare_test_cases()
