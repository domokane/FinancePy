# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.utils.date import Date

########################################################################################


def test_bump_bucket():

    value_dt = Date(1, 1, 2020)
    dates = [value_dt.add_years(1), value_dt.add_years(5), value_dt.add_years(10)]
    dfs = [0.97, 0.85, 0.70]
    bump = 0.0001

    curve = DiscountCurve(value_dt, dates, dfs)
    bumped = curve.bump_bucket(2.0, 4.0, bump)

    # Before the bucket the curve is unchanged
    assert np.isclose(bumped.df(dates[0]), 0.97, rtol=1e-12)

    # After the bucket the forward has been shifted only over [2, 4] years
    assert np.isclose(bumped.df(dates[1]), 0.85 * np.exp(-bump * 2.0), rtol=1e-12)
    assert np.isclose(bumped.df(dates[2]), 0.70 * np.exp(-bump * 2.0), rtol=1e-12)
