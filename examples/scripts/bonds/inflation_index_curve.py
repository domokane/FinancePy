# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import add_fp_to_path

from financepy.utils.date import Date
from financepy.products.inflation.inflation_index_curve import InflationIndexCurve



########################################################################################


def test_fin_inflation_index_curve():

    # Create a curve from times and discount factors
    index_dates = [Date(15, 1, 2008), Date(1, 4, 2008), Date(1, 5, 2008)]
    index_values = [209.49645, 214.823, 216.632]
    lag = 3  # months

    curve = InflationIndexCurve(index_dates, index_values, lag)

    ref_date = Date(22, 7, 2008)

    print("LABEL", "VALUE")

    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)

    print(ref_date, value)

    index_ratio = curve.index_ratio(ref_date)
    print(ref_date, index_ratio)


########################################################################################

#    print(curve)


test_fin_inflation_index_curve()
