########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import numpy as np
from scipy.stats import norm

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

value_dt = Date(1, 1, 2026)
expiry_dt = Date(1, 1, 2027)
t_exp = (expiry_dt - value_dt) / 365.0
r_dom, r_for, volatility, spot_fx_rate = 0.03, 0.05, 0.12, 1.10
lower_strike, upper_strike = 1.05, 1.20
domestic_curve = FlatDiscountCurve(value_dt, r_dom)
foreign_curve = FlatDiscountCurve(value_dt, r_for)
model = BlackScholes(volatility)

########################################################################################


def _double_digital(prem_currency, lower=lower_strike, upper=upper_strike):
    option = FXDoubleDigitalOption(expiry_dt, upper, lower, "EURUSD", 1.0, prem_currency)
    return option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)


def _digital_call(strike, prem_currency):
    option = FXDigitalOption(
        expiry_dt, strike, "EURUSD", OptionTypes.DIGITAL_CALL, 1.0, prem_currency
    )
    return option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)


def test_double_digital_is_a_difference_of_two_digital_calls():
    """Paying out between the two strikes is a digital call struck at the lower
    strike less one struck at the upper strike, in either premium currency."""
    for prem_currency in ("USD", "EUR"):
        expected = _digital_call(lower_strike, prem_currency) - _digital_call(
            upper_strike, prem_currency
        )
        assert np.isclose(_double_digital(prem_currency), expected, rtol=1e-10)


def test_double_digital_closed_form():
    sqrt_t = np.sqrt(t_exp)
    d2 = lambda k: (np.log(spot_fx_rate / k) + (r_dom - r_for - 0.5 * volatility**2) * t_exp) / (
        volatility * sqrt_t
    )
    d1 = lambda k: d2(k) + volatility * sqrt_t
    # domestic cash: discounted at the domestic rate, probabilities under d2
    domestic = np.exp(-r_dom * t_exp) * (norm.cdf(d2(lower_strike)) - norm.cdf(d2(upper_strike)))
    # one unit of foreign currency: asset-or-nothing, probabilities under d1
    foreign = (
        spot_fx_rate
        * np.exp(-r_for * t_exp)
        * (norm.cdf(d1(lower_strike)) - norm.cdf(d1(upper_strike)))
    )
    assert np.isclose(_double_digital("USD"), domestic, rtol=1e-6)
    assert np.isclose(_double_digital("EUR"), foreign, rtol=1e-6)


def test_double_digital_with_far_strikes_pays_with_certainty():
    """With the range covering every outcome the option is a certain payment: a
    domestic discount factor, or the spot value of one discounted foreign unit."""
    assert np.isclose(_double_digital("USD", 1e-6, 1e6), domestic_curve.df(expiry_dt))
    assert np.isclose(
        _double_digital("EUR", 1e-6, 1e6), spot_fx_rate * foreign_curve.df(expiry_dt)
    )
