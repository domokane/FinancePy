###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux
###############################################################################

import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.models.black import Black, BlackTypes, implied_volatility
from financepy.utils.date import Date
from financepy.utils.global_vars import gDaysInYear

forward = 0.034
strike = 0.050
riskFreeIR = 0.00
time_to_expiry = 2.0
volatility = 0.20
num_steps_per_year = 252.0

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


def test_implied_volatility():
    # European Exercise case
    # Input parametes are same as the Option on Treasury Futures Contract in the below link
    # http://gouthamanbalaraman.com/blog/value-options-commodity-futures-black-formula-quantlib-python.html
    strike = 119.0
    forward = 126.953
    time_to_expiry = 0.50
    volatility = 11.567/100.
    time_to_expiry = (Date(24, 12, 2015) - Date(1, 12, 2015)) / gDaysInYear
    r = 0.00105
    price = 7.968569878525531
    df = np.exp(-r*time_to_expiry)
    model = Black(volatility)
    # check value just in case
    valueCall = model.value(
        forward, strike, time_to_expiry, df, call_optionType)
    assert valueCall == price
    # check implied vol
    sigma = implied_volatility(
        forward, time_to_expiry, r, strike, price, OptionTypes.EUROPEAN_CALL)
    assert round(sigma, 5) == volatility
    # AmericanOptionImpliedVolatility: Implied Volatility calculation for American Option
    # https://rdrr.io/rforge/RQuantLib/man/AmericanOptionImpliedVolatility.html
    ...


def test_american_value_greeks():
    # Just check crr_tree_val_avg can called correctly from Black model
    strike = 100.0
    forward = 100.0
    time_to_expiry = 0.50
    volatility = 0.15
    modelTree = Black(volatility, implementation_type=BlackTypes.CRR_TREE,
                      num_steps=num_steps_per_year)
    expected = {
        'value': 4.229429096724559,
        'delta': -0.4788528545163757,
        'gamma': 0.037634509512885564,
        'theta': -4.233898040175347,
        'vega': 28.167013009550956
    }
    valueAmericanPut = modelTree.value(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT)
    assert valueAmericanPut == expected['value']
    deltaAmericanPut = modelTree.delta(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT)
    assert deltaAmericanPut == expected['delta']
    gammaAmericanPut = modelTree.gamma(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT)
    assert gammaAmericanPut == expected['gamma']
    thetaAmericanPut = modelTree.theta(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT)
    assert thetaAmericanPut == expected['theta']
    vegaAmericanPut = modelTree.vega(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT)
    assert vegaAmericanPut == expected['vega']
