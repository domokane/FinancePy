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


def test_default_floor_preserves_exact_expiry_for_scalar_and_list():
    """Exact zero time is used by existing deterministic expiry branches."""
    value_dt = Date(1, 1, 2026)
    assert option_years(value_dt, value_dt) == 0.0
    np.testing.assert_array_equal(option_years(value_dt, [value_dt]), [0.0])


def test_fail_false_preserves_negative_times_without_explicit_floor():
    """Opting out of validation does not implicitly clip expired dates."""
    value_dt = Date(1, 1, 2026)
    expiry_dt = Date(1, 1, 2025)
    assert option_years(value_dt, expiry_dt, fail=False) == -1.0
    np.testing.assert_array_equal(
        option_years(value_dt, [expiry_dt], fail=False), [-1.0]
    )


def test_index_option_at_expiry_keeps_zero_at_the_money_value():
    """Restoring vector support must not mask the merged index expiry fix."""
    from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
    from financepy.models.black import Black
    from financepy.products.equity.equity_index_option import EquityIndexOption
    from financepy.utils.global_types import OptionTypes

    value_dt = Date(1, 1, 2026)
    option = EquityIndexOption(value_dt, 100.0, OptionTypes.EUROPEAN_CALL)
    curve = FlatDiscountCurve(value_dt, 0.03)
    assert option.value(value_dt, 100.0, curve, Black(0.2)) == 0.0
