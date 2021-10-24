###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_zeros import DiscountCurveZeros
from financepy.market.curves.interpolator import InterpTypes
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
import numpy as np


def test_FinDiscountCurveZeros():
    start_date = Date(1, 1, 2018)
    times = np.linspace(1.0, 10.0, 10)
    dates = start_date.add_years(times)
    zero_rates = np.linspace(5.0, 6.0, 10)/100
    freq_type = FrequencyTypes.ANNUAL
    day_count_type = DayCountTypes.ACT_ACT_ISDA

    curve = DiscountCurveZeros(start_date,
                               dates,
                               zero_rates,
                               freq_type,
                               day_count_type,
                               InterpTypes.FLAT_FWD_RATES)

    date = start_date.add_years(0)
    df = curve.df(date)
    assert round(df, 4) == 1.0106

    date = start_date.add_years(2.5)
    df = curve.df(date)
    assert round(df, 4) == 0.8816

    date = start_date.add_years(5)
    df = curve.df(date)
    assert round(df, 4) == 0.7672

    date = start_date.add_years(7.5)
    df = curve.df(date)
    assert round(df, 4) == 0.6588

    date = start_date.add_years(10)
    df = curve.df(date)
    assert round(df, 4) == 0.5584
