###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
import sys
sys.path.append("..")


test_cases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinFlatCurve():

    curve_dt = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_dt.add_months(months)
    test_cases.header("COMPOUNDING", "DFS")
    compounding = FrequencyTypes.CONTINUOUS

    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    test_cases.print(compounding, dfs)

    compounding = FrequencyTypes.ANNUAL
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    test_cases.print(compounding, dfs)

    compounding = FrequencyTypes.SEMI_ANNUAL
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    test_cases.print(compounding, dfs)

    compounding = FrequencyTypes.QUARTERLY
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    test_cases.print(compounding, dfs)

    compounding = FrequencyTypes.MONTHLY
    flat_curve = DiscountCurveFlat(curve_dt, 0.05, compounding)
    dfs = flat_curve.df(dates)
    test_cases.print(compounding, dfs)

###############################################################################


test_FinFlatCurve()
test_cases.compareTestCases()
