###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date, set_date_format, DateFormatTypes
from financepy.products.rates.ibor_future import IborFuture

set_date_format(DateFormatTypes.UK_LONG)


def test_FinIborFuture():

    todayDate = Date(5, 5, 2020)

    i = 1
    fut = IborFuture(todayDate, i, "3M")
    fra = fut.to_fra(0.020, 0.0)
    assert fut._delivery_date == Date(17, 6, 2020)
    assert fra._start_date == Date(17, 6, 2020)

    i = 4
    fut = IborFuture(todayDate, i, "3M")
    fra = fut.to_fra(0.020, 0.0)
    assert fut._delivery_date == Date(17, 3, 2021)
    assert fra._start_date == Date(17, 3, 2021)

    i = 7
    fut = IborFuture(todayDate, i, "3M")
    fra = fut.to_fra(0.020, 0.0)
    assert fut._delivery_date == Date(15, 12, 2021)
    assert fra._start_date == Date(15, 12, 2021)

    i = 10
    fut = IborFuture(todayDate, i, "3M")
    fra = fut.to_fra(0.020, 0.0)
    assert fut._delivery_date == Date(21, 9, 2022)
    assert fra._start_date == Date(21, 9, 2022)
