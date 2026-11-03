# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
import pytest

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.models.heston import Heston, HestonValueTypes


# Reference see table 4.1 of Rouah book
value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 4, 2015)
v0 = 0.05  # initial variance of volatility
theta = 0.05  # long term variance
kappa = 2.0  # speed of variance reversion
sigma = 0.10  # volatility of variance
rho = -0.9  # correlation
interest_rate = 0.05
dividend_yield = 0.01
seed = 2838

num_steps = 100
num_paths = 20000
stock_price = 100.0

########################################################################################


def test_heston():

    rho = -0.90000
    sigma = 0.75000
    strike_price = 105.00

    heston_model = Heston(
        v0,
        kappa,
        theta,
        sigma,
        rho,
    )

    t_exp = (expiry_dt - value_dt) / 365

    opt_type = OptionTypes.EUROPEAN_CALL.value

    value_mc_heston = heston_model.value_mc(
        stock_price,
        t_exp,
        strike_price,
        opt_type,
        interest_rate,
        dividend_yield,
        num_paths,
        num_steps,
        seed,
    )

    value_gatheral = heston_model.value_call_gatheral(
        t_exp,
        strike_price,
        stock_price,
        interest_rate,
        dividend_yield,
    )

    value_lewis_rouah = heston_model.value_call_lewis_rouah(
        t_exp,
        strike_price,
        stock_price,
        interest_rate,
        dividend_yield,
    )

    value_lewis = heston_model.value_call_lewis(
        t_exp,
        strike_price,
        stock_price,
        interest_rate,
        dividend_yield,
    )

    value_weber = heston_model.value_call_weber(
        t_exp,
        strike_price,
        stock_price,
        interest_rate,
        dividend_yield,
    )

    assert value_mc_heston == pytest.approx(1.8626, abs=5.0e-3)
    assert value_gatheral == pytest.approx(1.8416, abs=5.0e-3)
    assert value_lewis_rouah == pytest.approx(1.8416, abs=5.0e-3)
    assert value_lewis == pytest.approx(1.8416, abs=5.0e-3)
    assert value_weber == pytest.approx(1.8416, abs=5.0e-3)


########################################################################################
# Issue #430: with a long maturity and kappa - rho * sigma < 0, the Gatheral and Weber
# P1 integrals have a 1/u tail and a boundary layer near u = 0 that quadrature on
# [0, inf) missed, so both methods returned a negative call price while LEWIS and
# LEWIS_ROUAH agreed with an independent reference (24.67708045).
########################################################################################


def test_heston_gatheral_weber_long_maturity_high_vol_of_vol():

    heston_model = Heston(v0=0.04, kappa=0.5, theta=0.04, xi=4.0, rho=0.5)

    s0 = 100.0
    t_exp = 30.0
    strike = 100.0
    r = 0.02
    q = 0.01

    reference = 24.67708045

    for method in HestonValueTypes:
        value = heston_model.call_value(s0, t_exp, strike, r, q, method=method)
        assert value == pytest.approx(reference, abs=1.0e-6), method


def test_heston_analytical_methods_agree():
    # Parameter sets that include mean-fleeing share-measure variance (kappa < rho * xi),
    # short and long maturities, and strikes on both sides of the forward.
    cases = [
        (0.04, 0.5, 0.04, 4.0, 0.5, 30.0, 60.0),
        (0.04, 0.5, 0.04, 4.0, 0.5, 10.0, 150.0),
        (0.09, 1.0, 0.09, 3.0, 0.7, 5.0, 100.0),
        (0.02, 3.0, 0.05, 0.3, -0.7, 1.0 / 365.0, 101.0),
        (0.05, 2.0, 0.05, 0.75, -0.9, 0.25, 105.0),
        (0.2, 0.2, 0.2, 5.0, 0.9, 50.0, 100.0),
    ]

    s0 = 100.0
    r = 0.03
    q = 0.01

    for v0, kappa, theta, xi, rho, t_exp, strike in cases:
        heston_model = Heston(v0=v0, kappa=kappa, theta=theta, xi=xi, rho=rho)
        values = [
            heston_model.call_value(s0, t_exp, strike, r, q, method=method)
            for method in HestonValueTypes
        ]
        lower = max(0.0, s0 * np.exp(-q * t_exp) - strike * np.exp(-r * t_exp))
        upper = s0 * np.exp(-q * t_exp)
        for value in values:
            assert lower - 1.0e-9 <= value <= upper + 1.0e-9
            assert value == pytest.approx(values[0], abs=1.0e-6)
