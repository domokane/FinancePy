"""Scalar Newton metadata must match the objective and derivative calls."""

import pytest

from financepy.utils.solver_1d import newton


@pytest.mark.parametrize("use_halley", [False, True])
def test_exact_initial_root_preserves_output_and_skips_derivatives(use_halley):
    """An exact initial root needs one objective call and no Newton step."""
    calls = []

    def objective(x, target):
        calls.append(x)
        return x - target

    def derivative(x, target):
        pytest.fail("No derivative should be evaluated at an exact root")

    root, info = newton(
        objective, 2.0, args=(2.0,), fprime=derivative,
        fprime2=derivative if use_halley else None, full_output=True,
    )
    assert root == 2.0
    assert info.root == root
    assert info.converged
    assert info.iterations == 0
    assert info.function_calls == len(calls) == 1


def test_exact_root_after_one_step_preserves_call_counts():
    """A linear solve reaches an exact root before a second Newton step."""
    calls = []

    def objective(x):
        calls.append(("value", x))
        return x - 2.0

    def derivative(x):
        calls.append(("derivative", x))
        return 1.0

    root, info = newton(objective, 0.0, fprime=derivative, full_output=True)
    assert root == 2.0
    assert info.converged
    assert info.iterations == 1
    assert info.function_calls == len(calls) == 3


def test_default_output_at_exact_root_is_still_scalar():
    """Metadata remains opt-in for initial and later exact-root exits."""
    assert newton(lambda x: x - 2.0, 2.0, fprime=lambda x: 1.0) == 2.0
    assert newton(lambda x: x - 2.0, 0.0, fprime=lambda x: 1.0) == 2.0
