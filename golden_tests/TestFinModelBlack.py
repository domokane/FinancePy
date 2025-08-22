# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux

import sys

sys.path.append("..")

import numpy as np
from financepy.models.black import Black
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_black():

    forward = 0.034
    strike = 0.050
    risk_free_ir = 0.00
    time_to_expiry = 2.0
    volatility = 0.20

    test_cases.header("ITEM", "CALL", "PUT")

    call_option_type = OptionTypes.EUROPEAN_CALL
    put_option_type = OptionTypes.EUROPEAN_PUT

    df = np.exp(-risk_free_ir * time_to_expiry)
    model = Black(volatility)

    dp = 12  # Precision

    try:

        value_call = model.value(forward, strike, time_to_expiry, df, call_option_type)
        value_put = model.value(forward, strike, time_to_expiry, df, put_option_type)

        assert round((value_call - value_put), dp) == round(
            df * (forward - strike), dp
        ), "The method called 'value()' doesn't comply with Call-Put parity"

        test_cases.print("VALUE", value_call, value_put)

        delta_call = model.delta(forward, strike, time_to_expiry, df, call_option_type)
        delta_put = model.delta(forward, strike, time_to_expiry, df, put_option_type)

        assert (
            round((1 / df) * (delta_call - delta_put), dp) == 1.0
        ), "The method called 'delta()' doesn't comply with Call-put parity"

        test_cases.print("DELTA", delta_call, delta_put)

        gamma_call = model.gamma(forward, strike, time_to_expiry, df, call_option_type)
        gamma_put = model.gamma(forward, strike, time_to_expiry, df, put_option_type)

        assert (
            round(gamma_call - gamma_put, dp) == 0.0
        ), "The method called 'gamma()' doesn't comply with Call-Put parity"

        test_cases.print("GAMMA", gamma_call, gamma_put)

        theta_call = model.theta(forward, strike, time_to_expiry, df, call_option_type)
        theta_put = model.theta(forward, strike, time_to_expiry, df, put_option_type)

        assert round((theta_call - theta_put), dp) == round(
            (risk_free_ir * time_to_expiry) * (forward - strike) * df, dp
        ), "The method called 'theta()' doesn't comply with Call-Put parity"

        test_cases.print("THETA", theta_call, theta_put)

        vega_call = model.vega(forward, strike, time_to_expiry, df, call_option_type)
        vega_put = model.vega(forward, strike, time_to_expiry, df, put_option_type)

        assert (
            round(vega_call - vega_put, dp) == 0.0
        ), "The method called 'vega()' doesn't comply with Call-Put parity"

        test_cases.print("VEGA", vega_call, vega_put)

    except AssertionError as err:
        raise err


########################################################################################

test_black()
test_cases.compare_test_cases()
