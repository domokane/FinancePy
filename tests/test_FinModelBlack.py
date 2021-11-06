###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux
###############################################################################

from financepy.utils.global_types import OptionTypes
from financepy.models.black import Black
import sys
import numpy as np

forward = 0.034
strike = 0.050
riskFreeIR = 0.00
time_to_expiry = 2.0
volatility = 0.20

call_optionType = OptionTypes.EUROPEAN_CALL
put_optionType = OptionTypes.EUROPEAN_PUT

df = np.exp(-riskFreeIR * time_to_expiry)
model = Black(volatility)

dp = 12  # Precision


def test_value():
    valueCall = model.value(
        forward, strike, time_to_expiry, df, call_optionType)
    valuePut = model.value(
        forward, strike, time_to_expiry, df, put_optionType)

    assert round((valueCall - valuePut), dp) == round(df*(forward - strike), dp), \
        "The method called 'value()' doesn't comply with Call-Put parity"

    assert round(valueCall*1000, 4) == 0.4599
    assert round(valuePut*10, 4) == 0.1646


def test_delta():
    deltaCall = model.delta(
        forward, strike, time_to_expiry, df, call_optionType)
    deltaPut = model.delta(
        forward, strike, time_to_expiry, df, put_optionType)

    assert round((1/df) * (deltaCall - deltaPut), dp) == 1.0, \
        "The method called 'delta()' doesn't comply with Call-put parity"

    assert round(deltaCall, 4) == 0.1108
    assert round(deltaPut, 4) == -0.8892


def test_gamma():
    gammaCall = model.gamma(
        forward, strike, time_to_expiry, df, call_optionType)
    gammaPut = model.gamma(
        forward, strike, time_to_expiry, df, put_optionType)

    assert round(gammaCall - gammaPut, dp) == 0.0, \
        "The method called 'gamma()' doesn't comply with Call-Put parity"

    assert round(gammaCall, 4) == 19.6594
    assert round(gammaPut, 4) == 19.6594


def test_theta():
    thetaCall = model.theta(
        forward, strike, time_to_expiry, df, call_optionType)
    thetaPut = model.theta(
        forward, strike, time_to_expiry, df, put_optionType)

    assert round((thetaCall - thetaPut), dp) == round((riskFreeIR * time_to_expiry) * (forward - strike) * df, dp), \
        "The method called 'theta()' doesn't comply with Call-Put parity"

    assert round(thetaCall*1000, 4) == -0.4545
    assert round(thetaPut*1000, 4) == -0.4545


def test_vega():
    vegaCall = model.vega(
        forward, strike, time_to_expiry, df, call_optionType)
    vegaPut = model.vega(
        forward, strike, time_to_expiry, df, put_optionType)

    assert round(vegaCall - vegaPut, dp) == 0.0, \
        "The method called 'vega()' doesn't comply with Call-Put parity"

    assert round(vegaCall*10, 4) == 0.0909
    assert round(vegaPut*10, 4) == 0.0909
