########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux
########################################################################################

import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.models.black import (
    Black,
    BlackTypes,
    black_value,
    implied_volatility,
)
from financepy.utils.date import Date
from financepy.utils.global_vars import G_DAYS_IN_YEARS
from financepy.models.equity_crr_tree import crr_tree_val_avg

forward = 0.034
strike = 0.050
risk_free_rate = 0.00
time_to_expiry = 2.0
volatility = 0.20
num_steps_per_year = 252.0

call_optionType = OptionTypes.EUROPEAN_CALL
put_optionType = OptionTypes.EUROPEAN_PUT

df = np.exp(-risk_free_rate * time_to_expiry)
model = Black(volatility)

dp = 12  # Precision


def test_value():
    value_call = model.value(
        forward, strike, time_to_expiry, df, call_optionType
    )
    value_put = model.value(
        forward, strike, time_to_expiry, df, put_optionType
    )

    assert round((value_call - value_put), dp) == round(
        df * (forward - strike), dp
    ), "The method called 'value()' doesn't comply with Call-Put parity"

    assert round(value_call * 1000, 4) == 0.4599
    assert round(value_put * 10, 4) == 0.1646


def test_delta():
    delta_call = model.delta(
        forward, strike, time_to_expiry, df, call_optionType
    )
    delta_put = model.delta(
        forward, strike, time_to_expiry, df, put_optionType
    )

    assert (
        round((1 / df) * (delta_call - delta_put), dp) == 1.0
    ), "The method called 'delta()' doesn't comply with Call-put parity"

    assert round(delta_call, 4) == 0.1108
    assert round(delta_put, 4) == -0.8892


def test_gamma():
    gamma_call = model.gamma(
        forward, strike, time_to_expiry, df, call_optionType
    )
    gamma_put = model.gamma(
        forward, strike, time_to_expiry, df, put_optionType
    )

    assert (
        round(gamma_call - gamma_put, dp) == 0.0
    ), "The method called 'gamma()' doesn't comply with Call-Put parity"

    assert round(gamma_call, 4) == 19.6594
    assert round(gamma_put, 4) == 19.6594


def test_theta():
    theta_call = model.theta(
        forward, strike, time_to_expiry, df, call_optionType
    )
    theta_put = model.theta(
        forward, strike, time_to_expiry, df, put_optionType
    )

    assert round((theta_call - theta_put), dp) == round(
        (risk_free_rate * time_to_expiry) * (forward - strike) * df, dp
    ), "The method called 'theta()' doesn't comply with Call-Put parity"

    assert round(theta_call * 1000, 4) == -0.4545
    assert round(theta_put * 1000, 4) == -0.4545


def test_vega():
    vega_call = model.vega(
        forward, strike, time_to_expiry, df, call_optionType
    )
    vega_put = model.vega(forward, strike, time_to_expiry, df, put_optionType)

    assert (
        round(vega_call - vega_put, dp) == 0.0
    ), "The method called 'vega()' doesn't comply with Call-Put parity"

    assert round(vega_call * 10, 4) == 0.0909
    assert round(vega_put * 10, 4) == 0.0909


def test_american_value_greeks():
    # Just check crr_tree_val_avg can called correctly from Black model
    # Expected result are obtained by derectly calling crr_tree_val_avg
    strike = 100.0
    forward = 100.0
    time_to_expiry = 0.50
    volatility = 0.15
    modelTree = Black(
        volatility,
        implementation_type=BlackTypes.CRR_TREE,
        num_steps=num_steps_per_year,
    )
    expected = {
        "value": 4.229429096724559,
        "delta": -0.4788528545163757,
        "gamma": 0.037634509512885564,
        "theta": -4.233898040175347,
        "vega": 28.167013009550956,
    }
    value_americanPut = modelTree.value(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT
    )
    assert round(value_americanPut, 5) == round(expected["value"], 5)
    deltaAmericanPut = modelTree.delta(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT
    )
    assert round(deltaAmericanPut, 5) == round(expected["delta"], 5)
    gammaAmericanPut = modelTree.gamma(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT
    )
    assert round(gammaAmericanPut, 5) == round(expected["gamma"], 5)
    thetaAmericanPut = modelTree.theta(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT
    )
    assert round(thetaAmericanPut, 5) == round(expected["theta"], 5)
    vegaAmericanPut = modelTree.vega(
        forward, strike, time_to_expiry, df, OptionTypes.AMERICAN_PUT
    )
    assert round(vegaAmericanPut, 5) == round(expected["vega"], 5)


def test_implied_volatility():
    # Implied Volatility calculation for European/American Call/Put Option
    # Just check if the method can reproduce the volatility of each model
    strike = 100.0
    forward = 100.0
    time_to_expiry = 0.50
    volatility = 0.15
    r = 0.1

    # European Call
    opt_type = OptionTypes.EUROPEAN_CALL
    price = black_value(
        forward, time_to_expiry, strike, r, volatility, opt_type
    )
    sigma = implied_volatility(
        forward, time_to_expiry, r, strike, price, opt_type
    )
    assert round(sigma, 5) == volatility

    # European Put
    opt_type = OptionTypes.EUROPEAN_PUT
    price = black_value(
        forward, time_to_expiry, strike, r, volatility, opt_type
    )
    sigma = implied_volatility(
        forward, time_to_expiry, r, strike, price, opt_type
    )
    assert round(sigma, 5) == volatility

    # American Call
    opt_type = OptionTypes.AMERICAN_CALL
    results = crr_tree_val_avg(
        forward,
        0.0,
        0.0,
        volatility,
        200,
        time_to_expiry,
        opt_type.value,
        strike,
    )
    price = results["value"]
    sigma = implied_volatility(
        forward, time_to_expiry, r, strike, price, opt_type
    )
    assert round(sigma, 5) == volatility

    # American Put
    opt_type = OptionTypes.AMERICAN_PUT
    results = crr_tree_val_avg(
        forward,
        0.0,
        0.0,
        volatility,
        200,
        time_to_expiry,
        opt_type.value,
        strike,
    )
    price = results["value"]
    sigma = implied_volatility(
        forward, time_to_expiry, r, strike, price, opt_type
    )
    assert round(sigma, 5) == volatility
