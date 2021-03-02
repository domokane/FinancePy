###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.utils.Date import Date
from financepy.utils.Frequency import FinFrequencyTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinFlatCurve():

    curveDate = Date(1, 1, 2019)
    months = range(1, 60, 3)
    dates = curveDate.addMonths(months)
    testCases.header("COMPOUNDING", "DFS")
    compounding = FinFrequencyTypes.CONTINUOUS

    flatCurve = FinDiscountCurveFlat(curveDate, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FinFrequencyTypes.ANNUAL
    flatCurve = FinDiscountCurveFlat(curveDate, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FinFrequencyTypes.SEMI_ANNUAL
    flatCurve = FinDiscountCurveFlat(curveDate, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FinFrequencyTypes.QUARTERLY
    flatCurve = FinDiscountCurveFlat(curveDate, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

    compounding = FinFrequencyTypes.MONTHLY
    flatCurve = FinDiscountCurveFlat(curveDate, 0.05, compounding)
    dfs = flatCurve.df(dates)
    testCases.print(compounding, dfs)

###############################################################################


test_FinFlatCurve()
testCases.compareTestCases()
