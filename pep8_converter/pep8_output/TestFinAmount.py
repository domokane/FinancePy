# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys

sys.path.append("..")

from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

################################################################################


def test_amount():

    test_cases.header("LABEL", "AMOUNT")
    x = Amount(101000.232)

    test_cases.print("Amount", x)

    x = Amount(101000.232, CurrencyTypes.CAD)

    test_cases.print("Amount", x)




################################################################################

test_amount()

test_cases.compare_test_cases()
