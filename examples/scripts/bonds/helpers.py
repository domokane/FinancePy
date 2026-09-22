#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import pandas as pd

# ============================================================================
# FINANCEPY EXAMPLES - Helpers
# ============================================================================


def parse_date_or_tenor(value):
    """Convert a DD/MM/YYYY value to Timestamp while preserving tenor strings.

    Benchmark data can contain either explicit dates or tenor strings such
    as "3M", "6M" and "5Y". Dates are converted to pandas Timestamp objects,
    while values that cannot be parsed as dates are returned unchanged.

    This replaces the deprecated pandas pattern:

        pd.to_datetime(..., errors="ignore")
    """

    if pd.isna(value):
        return value

    try:
        return pd.to_datetime(
            value,
            format="%d/%m/%Y",
        )
    except (ValueError, TypeError):
        return value
