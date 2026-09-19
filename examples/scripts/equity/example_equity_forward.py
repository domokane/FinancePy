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

from financepy.products.equity.equity_forward import EquityForward
from financepy.utils.date import Date
from financepy.utils.global_types import LongShortTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


########################################################################################


def test_equity_forward():

    value_dt = Date(13, 2, 2018)
    expiry_dt = value_dt.add_months(12)

    stock_price = 130.0
    forward_price = 125.0  # Locked
    discount_rate = 0.05
    dividend_rate = 0.02

    expiry_dt = value_dt.add_months(12)
    notional = 100.0

    discount_curve = FlatDiscountCurve(value_dt, discount_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

    equity_forward = EquityForward(expiry_dt, forward_price, notional, LongShortTypes.LONG)

    print("SPOT FX", "FX FWD", "VALUE_BS")

    fwd_price = equity_forward.forward(value_dt, stock_price, discount_curve, dividend_curve)

    fwd_value = equity_forward.value(value_dt, stock_price, discount_curve, dividend_curve)

    #    print(stock_price, fwd_price, fwd_value)
    print(stock_price, fwd_price, fwd_value)


########################################################################################

test_equity_forward()
