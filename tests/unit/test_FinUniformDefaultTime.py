import numpy as np
import pytest

from financepy.utils.helpers import uniform_to_default_time


def test_uniform_to_default_time_extrapolates_with_last_forward_hazard():
    times = np.array([0.0, 1.0, 2.0])
    survival = np.array([1.0, 0.9, 0.8])

    default_time = uniform_to_default_time(0.5, times, survival)

    expected = 2.0 + np.log(0.8 / 0.5) / np.log(0.9 / 0.8)
    assert default_time == pytest.approx(expected)


def test_uniform_to_default_time_interpolates_inside_curve():
    times = np.array([0.0, 1.0, 2.0])
    survival = np.array([1.0, 0.9, 0.8])

    default_time = uniform_to_default_time(0.85, times, survival)

    expected = 1.0 + np.log(0.9 / 0.85) / np.log(0.9 / 0.8)
    assert default_time == pytest.approx(expected)


def test_uniform_to_default_time_returns_no_default_for_flat_tail():
    times = np.array([0.0, 1.0, 2.0])
    survival = np.array([1.0, 0.9, 0.9])

    assert uniform_to_default_time(0.5, times, survival) == 99999.0
