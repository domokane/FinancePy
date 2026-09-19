# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path
_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause
_install_double_click_pause()
import numpy as np

import add_fp_to_path

from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date


########################################################################################


def test_fin_fx_digital_option():

    # Not exactly T=1.0 but close so don't exact exact agreement
    # (in fact I do not get exact agreement even if I do set T=1.0)
    value_dt = Date(13, 2, 2018)
    expiry_dt = Date(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1_cc_rate = 0.030  # EUR
    ccy2_cc_rate = 0.025  # USD

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 1.20
    strike_fx_rate = 1.250
    volatility = 0.10

    notional = 1.0

    domestic_curve = FlatDiscountCurve(value_dt, ccy2_cc_rate)
    foreign_curve = FlatDiscountCurve(value_dt, ccy1_cc_rate)

    model = BlackScholes(volatility)

    digital_option = FXDigitalOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.DIGITAL_CALL,
        notional,
        "USD",
    )

    spot_fx_rate = np.linspace(0.01, 2.0, 10)

    value = digital_option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)


########################################################################################

test_fin_fx_digital_option()
