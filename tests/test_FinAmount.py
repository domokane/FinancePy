##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount

##########################################################################


def test_Amount():

    amt = Amount(101000.232)
    assert amt.amount() == 101000.232

    amt = Amount(101000.232, CurrencyTypes.CAD)
    assert repr(amt) == "CAD 101,000.23"


###############################################################################
