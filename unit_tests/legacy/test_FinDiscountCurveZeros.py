# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes

########################################################################################


def test_fin_discount_curve_zeros():

    start_dt = Date(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    dates = start_dt.add_years(times)
    zero_rates = np.linspace(5.0, 6.0, 10) / 100
    freq_type = FrequencyTypes.ANNUAL

    curve = DiscountCurveZeros(
        start_dt,
        dates,
        zero_rates,
        freq_type,
        InterpTypes.FLAT_FWD_RATES,
    )

    date = start_dt.add_years(0)
    df = curve.df(date)
    assert round(df, 4) == 1.0106

    date = start_dt.add_years(2.5)
    df = curve.df(date)
    assert round(df, 4) == 0.8815

    date = start_dt.add_years(5)
    df = curve.df(date)
    assert round(df, 4) == 0.767

    date = start_dt.add_years(7.5)
    df = curve.df(date)
    assert round(df, 4) == 0.6586

    date = start_dt.add_years(10)
    df = curve.df(date)
    assert round(df, 4) == 0.5582
