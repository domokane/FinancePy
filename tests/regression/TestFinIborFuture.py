# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.products.rates.ibor_future import IborFuture
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

set_date_format(DateFormatTypes.UK_LONG)

########################################################################################


def test_fin_ibor_future():

    today_date = Date(5, 5, 2020)

    test_cases.header("VALUES")

    for i in range(1, 12):
        fut = IborFuture(today_date, i, "3M")
        test_cases.print(fut)

        fra = fut.to_fra(0.020, 0.0)
        test_cases.print(fra)


########################################################################################

test_fin_ibor_future()
test_cases.compare_test_cases()
