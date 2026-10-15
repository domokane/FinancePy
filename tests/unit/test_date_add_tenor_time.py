"""Regressions for intraday Date.add_tenor consistency (issue #390)."""

import pytest

from financepy.utils.date import Date


@pytest.mark.parametrize(
    ("tenor", "ymd"),
    [
        ("1D", (3, 1, 2026)),
        ("-1D", (1, 1, 2026)),
        ("0D", (2, 1, 2026)),
        ("1W", (9, 1, 2026)),
        ("-1W", (26, 12, 2025)),
        ("0W", (2, 1, 2026)),
        ("1M", (2, 2, 2026)),
        ("-1M", (2, 12, 2025)),
        ("0M", (2, 1, 2026)),
        ("1Y", (2, 1, 2027)),
        ("-1Y", (2, 1, 2025)),
        ("0Y", (2, 1, 2026)),
    ],
)
def test_add_tenor_preserves_intraday_fields_and_serial(tenor, ymd):
    original = Date(2, 1, 2026, 12, 34, 56)
    result = original.add_tenor(tenor)
    expected = Date(*ymd, 12, 34, 56)

    assert result is not original
    assert result.datetime() == expected.datetime()
    assert result.excel_dt == pytest.approx(expected.excel_dt, rel=0, abs=1e-10)
    assert original.datetime() == Date(2, 1, 2026, 12, 34, 56).datetime()


def test_add_tenor_list_retains_time_for_each_unit():
    original = Date(2, 1, 2026, 12, 34, 56)
    results = original.add_tenor(["1D", "1W", "1M", "1Y"])

    assert [result.datetime() for result in results] == [
        Date(3, 1, 2026, 12, 34, 56).datetime(),
        Date(9, 1, 2026, 12, 34, 56).datetime(),
        Date(2, 2, 2026, 12, 34, 56).datetime(),
        Date(2, 1, 2027, 12, 34, 56).datetime(),
    ]


def test_add_tenor_month_end_preserves_time_and_leap_day_clamping():
    month_end = Date(31, 1, 2024, 23, 59, 59).add_tenor("1M")
    leap_day = Date(29, 2, 2024, 9, 8, 7).add_tenor("1Y")

    assert month_end.datetime() == Date(29, 2, 2024, 23, 59, 59).datetime()
    assert leap_day.datetime() == Date(28, 2, 2025, 9, 8, 7).datetime()
    assert month_end.excel_dt == pytest.approx(
        Date(29, 2, 2024, 23, 59, 59).excel_dt, rel=0, abs=1e-10
    )
