# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

import time
import numpy as np

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes


test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_discount_curve_zeros():

    start_dt = Date(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    dates = start_dt.add_years(times)
    zero_rates = np.linspace(5.0, 6.0, 10) / 100
    freq_type = FrequencyTypes.ANNUAL
    dc_type = DayCountTypes.ACT_ACT_ISDA

    curve = DiscountCurveZeros(
        start_dt,
        dates,
        zero_rates,
        freq_type,
        dc_type,
        InterpTypes.FLAT_FWD_RATES,
    )

    test_cases.header("T", "DF")

    years = np.linspace(0, 10, 21)
    dates = start_dt.add_years(years)
    for dt in dates:
        df = curve.df(dt)
        test_cases.print(dt, df)

    #    print(curve)

    num_repeats = 100

    start = time.time()

    for i in range(0, num_repeats):
        freq_type = FrequencyTypes.ANNUAL
        dc_type = DayCountTypes.ACT_ACT_ISDA

        dates = [
            Date(14, 6, 2016),
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
            0.000000,
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

        start_dt = dates[0]

        curve = DiscountCurveZeros(
            start_dt,
            dates,
            zero_rates,
            freq_type,
            dc_type,
            InterpTypes.FLAT_FWD_RATES,
        )

    end = time.time()
    period = end - start


########################################################################################

#    print("Time taken:", period)

#    print(curve)


test_fin_discount_curve_zeros()
test_cases.compare_test_cases()
