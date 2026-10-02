"""Nelder-Mead convergence requires a small simplex and objective spread."""

import numpy as np
import pytest
from numba import njit

from financepy.utils.solver_nm import nelder_mead


@njit
def quadratic(x, target, scale):
    """Quadratic with a known minimizer and independently variable scale."""
    return scale * np.sum((x - target) ** 2)


@pytest.mark.parametrize("target,scale", [(0.205, 1.0), (0.7, 1.0), (0.7, 1e-12)])
def test_equal_or_small_objective_spread_does_not_stop_wide_simplex(target, scale):
    """Equal endpoint values do not establish proximity to a minimum."""
    tol_x = 1e-8
    result = nelder_mead(
        quadratic, np.array([0.2]), args=(target, scale),
        tol_f=1e-10, tol_x=tol_x,
    )
    assert result.success
    assert abs(result.x[0] - target) < 2 * tol_x
    assert np.max(np.abs(result.final_simplex - result.x)) < tol_x


def test_large_coordinate_simplex_is_measured_in_coordinate_units():
    """A dimensionless volume ratio must not certify coordinate tolerance."""
    result = nelder_mead(
        quadratic, np.array([1e6, -1e6]), args=(700000.0, 1.0),
        tol_x=1e-6, tol_f=1e-8,
    )
    assert result.success
    assert np.max(np.abs(result.x - 700000.0)) < 2e-6
    assert np.max(np.abs(result.final_simplex - result.x)) < 1e-6


def test_iteration_limit_does_not_report_convergence():
    """An exhausted budget with a wide simplex remains unsuccessful."""
    result = nelder_mead(
        quadratic, np.array([0.2]), args=(0.7, 1.0), max_iter=0,
    )
    assert not result.success
    assert result.nit == 0
