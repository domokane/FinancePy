# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount

########################################################################################


def test__amount():

    amt = Amount(101000.232)
    assert amt.amount() == 101000.232

    amt = Amount(101000.232, CurrencyTypes.CAD)
    assert repr(amt) == "CAD 101,000.23"

