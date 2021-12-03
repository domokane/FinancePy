###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################


def test_FinInflationIndexCurve():

    # Create a curve from times and discount factors
    indexDates = [Date(15, 1, 2008), Date(1, 4, 2008), Date(1, 5, 2008)]
    indexValues = [209.49645, 214.823, 216.632]
    lag = 3  # months

    curve = FinInflationIndexCurve(indexDates, indexValues, lag)

    refDate = Date(22, 7, 2008)

    testCases.header("LABEL", "VALUE")

    value = curve.index_value(refDate)
    value = curve.index_value(refDate)
    value = curve.index_value(refDate)
    value = curve.index_value(refDate)

    testCases.print(refDate, value)

    index_ratio = curve.index_ratio(refDate)
    testCases.print(refDate, index_ratio)

#    print(curve)

###############################################################################


test_FinInflationIndexCurve()
testCases.compareTestCases()
