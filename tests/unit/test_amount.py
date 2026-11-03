# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.amount import Amount

########################################################################################


def test_amount():

    amt = Amount(101000.232, "CAD")
    assert repr(amt) == "CAD 101,000.23"
