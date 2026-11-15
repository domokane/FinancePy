# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.frequency import FrequencyTypes

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


@pytest.mark.parametrize("container", [list, tuple, np.asarray])
@pytest.mark.parametrize("frequency", [
    FrequencyTypes.CONTINUOUS,
    FrequencyTypes.SIMPLE_INTEREST,
    FrequencyTypes.ANNUAL,
])
def test_zero_rate_time_containers(container, frequency):
    curve = FlatDiscountCurve(Date(11, 11, 2026), 0.03)
    times = np.array([1.0, 2.0, 5.0])
    actual = curve.zero_rate_t(container(times), frequency)
    expected = curve.zero_rate_t(times, frequency)
    np.testing.assert_allclose(actual, expected, rtol=1e-12)
    assert isinstance(actual, np.ndarray)


@pytest.mark.parametrize("container", [list, tuple, np.asarray])
def test_zero_rate_negative_time_containers(container):
    curve = FlatDiscountCurve(Date(11, 11, 2026), 0.03)
    with pytest.raises(FinError, match="Times to maturity cannot be negative"):
        curve.zero_rate_cc_t(container([1.0, -1.0]))


def test_zero_rate_scalar_and_empty_inputs():
    curve = FlatDiscountCurve(Date(11, 11, 2026), 0.03)
    assert isinstance(curve.zero_rate_cc_t(1.0), float)
    assert curve.zero_rate_cc_t(1.0) == pytest.approx(0.03)
    np.testing.assert_array_equal(curve.zero_rate_cc_t([]), np.array([]))
