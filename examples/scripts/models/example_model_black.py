# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux


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

from financepy.models.black import Black
from financepy.utils.global_types import OptionTypes


########################################################################################


def test_black():

    forward = 0.034
    strike = 0.050
    risk_free_ir = 0.00
    t_exp = 2.0
    volatility = 0.20

    print("ITEM", "CALL", "PUT")

    call_option_type = OptionTypes.EUROPEAN_CALL
    put_option_type = OptionTypes.EUROPEAN_PUT

    df = np.exp(-risk_free_ir * t_exp)
    model = Black(volatility)

    dp = 12  # Precision

    try:

        value_call = model.value(forward, strike, t_exp, df, call_option_type)
        value_put = model.value(forward, strike, t_exp, df, put_option_type)

        assert round((value_call - value_put), dp) == round(
            df * (forward - strike), dp
        ), "The method called 'value()' doesn't comply with Call-Put parity"

        print("VALUE", value_call, value_put)

        delta_call = model.delta(forward, strike, t_exp, df, call_option_type)
        delta_put = model.delta(forward, strike, t_exp, df, put_option_type)

        assert (
            round((1 / df) * (delta_call - delta_put), dp) == 1.0
        ), "The method called 'delta()' doesn't comply with Call-put parity"

        print("DELTA", delta_call, delta_put)

        gamma_call = model.gamma(forward, strike, t_exp, df, call_option_type)
        gamma_put = model.gamma(forward, strike, t_exp, df, put_option_type)

        assert (
            round(gamma_call - gamma_put, dp) == 0.0
        ), "The method called 'gamma()' doesn't comply with Call-Put parity"

        print("GAMMA", gamma_call, gamma_put)

        theta_call = model.theta(forward, strike, t_exp, df, call_option_type)
        theta_put = model.theta(forward, strike, t_exp, df, put_option_type)

        assert round((theta_call - theta_put), dp) == round(
            (risk_free_ir * t_exp) * (forward - strike) * df, dp
        ), "The method called 'theta()' doesn't comply with Call-Put parity"

        print("THETA", theta_call, theta_put)

        vega_call = model.vega(forward, strike, t_exp, df, call_option_type)
        vega_put = model.vega(forward, strike, t_exp, df, put_option_type)

        assert (
            round(vega_call - vega_put, dp) == 0.0
        ), "The method called 'vega()' doesn't comply with Call-Put parity"

        print("VEGA", vega_call, vega_put)

    except AssertionError as err:
        raise err


########################################################################################

test_black()
