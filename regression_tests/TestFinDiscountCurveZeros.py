# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time
import numpy as np

import add_fp_to_path

from financepy.market.curves.zero_rates_discount_curve import ZeroRatesDiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes


from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_discount_curve_zeros():

    dates = [
        Date(14, 9, 2016),
        Date(14, 12, 2016),
        Date(14, 6, 2017),
        Date(14, 6, 2019),
        Date(14, 6, 2021),
        Date(15, 6, 2026),
        Date(16, 6, 2031),
        Date(16, 6, 2036),
        Date(14, 6, 2046),
    ]

    zero_rates = [
        0.006616,
        0.007049,
        0.007795,
        0.009599,
        0.011203,
        0.015068,
        0.017583,
        0.018998,
        0.020080,
    ]

    start_dt = Date(14, 6, 2016)

    times = np.linspace(0.0, 30, 100)

    test_cases.header("Interp_type", "Time", "Zero_cc", "Fwd_cc", "Calc_Time")

    for interp_type in InterpTypes:

        start = time.time()

        freq_type = FrequencyTypes.ANNUAL
        time_dc_type = DayCountTypes.ACT_ACT_ISDA

        curve = ZeroRatesDiscountCurve(
            start_dt,
            dates,
            zero_rates,
            freq_type,
            interp_type,
            time_dc_type,
        )

        zeros_cc = curve.zero_rate_cc_t(times) * 100.0
        fwd_cc = curve.fwd_rate_inst_t(times) * 100.0

        end = time.time()
        period = end - start

        for t, z, f in zip(times, zeros_cc, fwd_cc):
            test_cases.print(interp_type, t, z, f, period)


########################################################################################


test_fin_discount_curve_zeros()
test_cases.compare_test_cases()
