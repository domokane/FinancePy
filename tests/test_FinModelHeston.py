###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.utils.global_types import OptionTypes
from financepy.models.heston import Heston, HestonNumericalScheme
import numpy as np


# Reference see table 4.1 of Rouah book
valuation_date = Date(1, 1, 2015)
expiry_date = Date(1, 4, 2015)
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


def test_heston():
    rho = -0.90000
    sigma = 0.75000
    strike_price = 105.00
    hestonModel = Heston(v0, kappa, theta, sigma, rho)

    call_option = EquityVanillaOption(
        expiry_date, strike_price, OptionTypes.EUROPEAN_CALL)

    value_mc_Heston = hestonModel.value_mc(
        valuation_date,
        call_option,
        stock_price,
        interest_rate,
        dividend_yield,
        num_paths,
        num_steps,
        seed)
    valueGatheral = hestonModel.value_gatheral(
        valuation_date, call_option, stock_price, interest_rate, dividend_yield)
    valueLewisRouah = hestonModel.value_lewis_rouah(
        valuation_date, call_option, stock_price, interest_rate, dividend_yield)
    valueLewis = hestonModel.value_lewis(
        valuation_date, call_option, stock_price, interest_rate, dividend_yield)
    valueWeber = hestonModel.value_weber(
        valuation_date, call_option, stock_price, interest_rate, dividend_yield)

    assert round(value_mc_Heston, 4) == 1.7333
    assert round(valueGatheral, 4) == 1.8416
    assert round(valueLewisRouah, 4) == 1.8416
    assert round(valueLewis, 4) == 1.8416
    assert round(valueWeber, 4) == 1.8416
