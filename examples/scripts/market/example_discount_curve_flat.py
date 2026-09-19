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

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve




########################################################################################


def test_fin_flat_curve():

    curve_dt = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_dt.add_months(months)
    print("COMPOUNDING", "DFS")
    compounding = FrequencyTypes.CONTINUOUS

    flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    print(compounding, dfs)

    compounding = FrequencyTypes.ANNUAL
    flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    print(compounding, dfs)

    compounding = FrequencyTypes.SEMI_ANNUAL
    flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    print(compounding, dfs)

    compounding = FrequencyTypes.QUARTERLY
    flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    print(compounding, dfs)

    compounding = FrequencyTypes.MONTHLY
    flat_curve = FlatDiscountCurve(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    print(compounding, dfs)


########################################################################################

test_fin_flat_curve()
