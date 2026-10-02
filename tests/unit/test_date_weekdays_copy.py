import pytest

from financepy.utils.date import Date


@pytest.mark.parametrize("day", [2, 3, 4])
def test_zero_weekdays_preserves_weekday_and_weekend_dates_without_aliasing(day):
    original = Date(day, 1, 2026)
    result = original.add_weekdays(0)
    assert result == original
    assert result is not original
    result.d = 20
    assert original.d == day


def test_zero_weekdays_preserves_intraday_value():
    original = Date(2, 1, 2026, 12, 34, 56)
    result = original.add_weekdays(0)
    assert result is not original
    assert result.excel_dt == original.excel_dt
    assert (result.hh, result.mm, result.ss) == (12, 34, 56)
