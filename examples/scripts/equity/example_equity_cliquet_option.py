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

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes



########################################################################################


def test_equity_cliquet_option():

    start_dt = Date(1, 1, 2014)
    final_expiry_dt = Date(1, 1, 2017)
    freq_type = FrequencyTypes.QUARTERLY
    opt_type = OptionTypes.EUROPEAN_CALL

    cliquet_option = EquityCliquetOption(start_dt, final_expiry_dt, opt_type, freq_type)

    value_dt = Date(1, 1, 2015)
    stock_price = 100.0
    volatility = 0.20
    interest_rate = 0.05
    dividend_yield = 0.02
    model = BlackScholes(volatility)
    discount_curve = FlatDiscountCurve(value_dt, interest_rate)
    dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

    v = cliquet_option.value(
        value_dt, stock_price, discount_curve, dividend_curve, model
    )

    print("LABEL", "VALUE")
    print("FINANCEPY", v)


########################################################################################

test_equity_cliquet_option()
