# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.amount import Amount

########################################################################################


def test_amount():

    amt = Amount(101000.232, "CAD")
    assert repr(amt) == "CAD 101,000.23"


def test_amount_sub_and_eq():

    diff = Amount(100.0, "CAD") - Amount(40.0, "CAD")
    assert diff.amount == 60.0
    assert diff.ccy_name == "CAD"
    assert Amount(1.0, "CAD") == Amount(1.0, "CAD")
    assert Amount(1.0, "CAD") != Amount(1.0, "USD")
