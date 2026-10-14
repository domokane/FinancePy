"""Regression coverage for bounded Nelder-Mead objective evaluations."""

import numpy as np
import pytest
from numba import njit

from financepy.utils.solver_nm import nelder_mead


@njit
def bounded_quadratic(x, target, lower, upper):
    """Reject out-of-domain evaluations and minimize distance to the target."""
    if np.any(x < lower) or np.any(x > upper):
        raise ValueError("Objective evaluated outside bounds")
    return np.sum((x - target) ** 2)


@pytest.mark.parametrize(
    "initial,target,expected",
    [(1.0, 2.0, 1.0), (0.0, -1.0, 0.0), (3.0, 2.0, 1.0)],
)
def test_bounded_solver_never_evaluates_outside_domain(initial, target, expected):
    """Cover upper/lower optima and an out-of-domain initial guess."""
    bounds = np.array([[0.0, 1.0]])
    result = nelder_mead(
        bounded_quadratic,
        np.array([initial]),
        bounds=bounds,
        args=(np.array([target]), bounds[:, 0], bounds[:, 1]),
    )

    assert result.success
    assert result.x[0] == pytest.approx(expected, abs=1e-4)
    assert np.all(result.final_simplex >= bounds[:, 0])
    assert np.all(result.final_simplex <= bounds[:, 1])


def test_bounded_solver_handles_multiple_variables_and_function_arguments():
    """Check the constrained optimum of a two-variable objective."""
    bounds = np.array([[0.0, 1.0], [-1.0, 0.0]])
    target = np.array([2.0, -2.0])
    result = nelder_mead(
        bounded_quadratic,
        np.array([0.5, -0.5]),
        bounds=bounds,
        args=(target, bounds[:, 0], bounds[:, 1]),
    )

    assert result.success
    np.testing.assert_allclose(result.x, [1.0, -1.0], atol=1e-4)
    assert result.fun == pytest.approx(2.0, abs=1e-4)
