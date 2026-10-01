import pytest

from financepy.utils.date import Date
from financepy.utils.date_format import (
    DateFormatTypes,
    get_date_format,
    set_date_format,
)
from financepy.utils.error import FinError


def test_set_date_format_accepts_enum_name_and_keeps_enum_state():
    try:
        set_date_format("us_short")
        assert get_date_format() is DateFormatTypes.US_SHORT
        assert str(Date(12, 3, 2026)) == "03-12-26"

        set_date_format(DateFormatTypes.UK_LONG)
        assert str(Date(12, 3, 2026)) == "12-MAR-2026"
    finally:
        set_date_format(DateFormatTypes.UK_LONG)


def test_set_date_format_rejects_unknown_name():
    with pytest.raises(FinError, match="Unknown date format"):
        set_date_format("NOT_A_DATE_FORMAT")
