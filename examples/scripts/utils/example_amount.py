# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import add_fp_to_path

from financepy.utils.currency import CurrencyTypes
from financepy.utils.amount import Amount



########################################################################################


def test_amount():

    print("LABEL", "AMOUNT")
    x = Amount(101000.232, CurrencyTypes.USD)

    print("Amount", x)

    x = Amount(101000.232, CurrencyTypes.CAD)

    print("Amount", x)


########################################################################################

test_amount()

