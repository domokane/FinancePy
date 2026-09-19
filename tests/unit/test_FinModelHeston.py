# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from financepy.utils.date import Date
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.global_types import OptionTypes, HestonNumericalSchemeTypes
from financepy.models.heston import Heston

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


import pytest


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
