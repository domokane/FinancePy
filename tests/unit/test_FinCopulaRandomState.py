"""Regression tests for copula simulation random state isolation."""

from types import SimpleNamespace

import numpy as np
import pytest

from financepy.models.gauss_copula import default_times_gc
from financepy.models.student_t_copula import StudentTCopula


@pytest.mark.parametrize("model", ["gaussian", "student_t"])
def test_default_time_sampling_does_not_change_global_random_state(model):
    curve = SimpleNamespace(
        _times=np.array([0.0, 1.0, 2.0]),
        _qs=np.array([1.0, 0.8, 0.6]),
    )
    original_state = np.random.get_state()
    try:
        np.random.seed(19)
        expected_next = np.random.random()
        np.random.seed(19)
        if model == "gaussian":
            result = default_times_gc([curve], np.eye(1), num_trials=3, seed=42)
            expected = [
                4.298859600830989,
                2.2710975345549813,
                4.92570788159172,
                1.512644582129489,
                3.0387161160445206,
                1.2643966898636527,
            ]
        else:
            result = StudentTCopula().default_times(
                [curve], np.eye(1), degrees_of_freedom=5, num_trials=3, seed=42
            )
            expected = [
                1.8475215516073522,
                3.009593559034871,
                1.542702875015536,
                3.6514637162971484,
                2.2946270938944044,
                4.232781171714952,
            ]
        assert result.shape == (1, 6)
        np.testing.assert_allclose(result[0], expected, rtol=1e-10)
        assert np.random.random() == expected_next
    finally:
        np.random.set_state(original_state)
