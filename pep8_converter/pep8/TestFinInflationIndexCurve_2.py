# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

from financepy.utils.date import Date
from financepy.products.inflation.FinInflationIndexCurve import (
    FinInflationIndexCurve,
)
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)



########################################################################################


def test_FinInflationIndexCurve():

    # Create a curve from times and discount factors
    index_dates = [Date(15, 1, 2008), Date(1, 4, 2008), Date(1, 5, 2008)]
    index_values = [209.49645, 214.823, 216.632]
    lag = 3  # months

    curve = FinInflationIndexCurve(index_dates, index_values, lag)

    ref_date = Date(22, 7, 2008)

    test_cases.header("LABEL", "VALUE")

    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)
    value = curve.index_value(ref_date)

    test_cases.print(ref_date, value)

    index_ratio = curve.index_ratio(ref_date)
    test_cases.print(ref_date, index_ratio)


#    print(curve)

########################################################################################


test_FinInflationIndexCurve()
test_cases.compare_test_cases()

########################################################################################

