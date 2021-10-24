###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.products.inflation.FinInflationIndexCurve import FinInflationIndexCurve
from financepy.utils.date import Date


def test_FinInflationIndexCurve():
    # Create a curve from times and discount factors
    indexDates = [Date(15, 1, 2008), Date(1, 4, 2008), Date(1, 5, 2008)]
    indexValues = [209.49645, 214.823, 216.632]
    lag = 3  # months

    curve = FinInflationIndexCurve(indexDates, indexValues, lag)
    refDate = Date(22, 7, 2008)

    value = curve.index_value(refDate)
    assert round(value, 4) == 216.0485

    index_ratio = curve.index_ratio(refDate)
    assert round(index_ratio, 4) == 1.0313
