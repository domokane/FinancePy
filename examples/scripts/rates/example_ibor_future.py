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

from financepy.products.rates.ibor_future import IborFuture
from financepy.utils.date_format import set_date_format, DateFormatTypes
from financepy.utils.date import Date


set_date_format(DateFormatTypes.UK_LONG)

########################################################################################


def test_fin_ibor_future():

    today_date = Date(5, 5, 2020)

    print("VALUES")

    for i in range(1, 12):
        fut = IborFuture(today_date, i, "3M")
        print(fut)

        fra = fut.to_fra(0.020, 0.0)
        print(fra)


########################################################################################

test_fin_ibor_future()
