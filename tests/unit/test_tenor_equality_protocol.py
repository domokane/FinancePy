"""Tenors participate in ordinary heterogeneous Python comparisons."""

import pytest

from financepy.utils.tenor import Tenor


@pytest.mark.parametrize("other", [None, "1Y", 1, object()])
def test_unrelated_values_are_not_equal(other):
    tenor = Tenor("1Y")
    assert tenor != other
    assert other != tenor
    assert tenor not in [other]


def test_other_operand_can_handle_equality():
    class MatchesTenor:
        def __eq__(self, other):
            return isinstance(other, Tenor)

    assert Tenor("1Y") == MatchesTenor()
