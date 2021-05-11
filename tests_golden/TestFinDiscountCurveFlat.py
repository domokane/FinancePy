###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.market.curves.curve_flat import DiscountCurveFlat
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinFlatCurve():

    curve_date = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_date.add_months(months)
    testCases.header("COMPOUNDING", "DFS")
    compounding = FrequencyTypes.CONTINUOUS

    flatCurve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FrequencyTypes.ANNUAL
    flatCurve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FrequencyTypes.SEMI_ANNUAL
    flatCurve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FrequencyTypes.QUARTERLY
    flatCurve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FrequencyTypes.MONTHLY
    flatCurve = DiscountCurveFlat(curve_date, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

###############################################################################


test_FinFlatCurve()
testCases.compareTestCases()
