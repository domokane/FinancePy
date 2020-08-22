###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode


from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinFrequency import FinFrequencyTypes
import sys
sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)


def test_FinFlatCurve():

    curveDate = FinDate(1, 1, 2019)
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
