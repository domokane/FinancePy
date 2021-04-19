###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.market.curves.curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinFlatCurve():

    curve_date = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curve_date.addMonths(months)
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
