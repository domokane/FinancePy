import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.utils.error import FinError
from financepy.utils.helpers import option_years


def test_option_years_supports_expiry_date_lists_and_floor():
    value_dt = Date(1, 1, 2026)
    expiries = [Date(1, 1, 2025), value_dt, Date(1, 1, 2027)]

    result = option_years(value_dt, expiries, floor=0.1, fail=False)

    np.testing.assert_allclose(result, [0.1, 0.1, 1.0])


def test_option_years_rejects_expiry_before_value_date_in_list():
    value_dt = Date(1, 1, 2026)

    with pytest.raises(FinError, match="Option expires before value date"):
        option_years(value_dt, [Date(1, 1, 2027), Date(1, 1, 2025)])
