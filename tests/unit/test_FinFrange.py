import pytest

from financepy.utils.error import FinError
from financepy.utils.helpers import frange as helper_frange
from financepy.utils.math import frange as math_frange


@pytest.mark.parametrize("frange", [helper_frange, math_frange])
def test_frange_supports_steps_toward_stop(frange):
    assert frange(1, 5, 2) == [1, 3, 5]
    assert frange(5, 1, -2) == [5, 3, 1]


@pytest.mark.parametrize("frange", [helper_frange, math_frange])
def test_frange_rejects_zero_step(frange):
    with pytest.raises(FinError, match="Step cannot be zero"):
        frange(1, 5, 0)
