###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

from financepy.finutils.FinDate import FinDate
from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

##############################################################################

def test_FinInflationIndexCurve():

    # Create a curve from times and discount factors
    indexDates = [FinDate(15,1,2008), FinDate(1,4,2008), FinDate(1,5,2008)]
    indexValues = [209.49645, 214.823, 216.632]
    lag = 3 # months

    curve = FinInflationIndexCurve(indexDates, indexValues, lag)

    refDate = FinDate(22, 7, 2008)
    
    testCases.header("LABEL", "VALUE")

    value = curve.indexValue(refDate)
    value = curve.indexValue(refDate)
    value = curve.indexValue(refDate)
    value = curve.indexValue(refDate)

    testCases.print(refDate, value)

    indexRatio = curve.indexRatio(refDate)
    testCases.print(refDate, indexRatio)

#    print(curve)

###############################################################################


test_FinInflationIndexCurve()
testCases.compareTestCases()
