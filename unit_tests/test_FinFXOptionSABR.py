###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
import numpy as np


# UNFINISHED
# There is no FXAmericanOption class. It is embedded in the FXVanillaOption
# class. This test just compares it to the European

value_dt = Date(13, 2, 2018)
expiry_dt = Date(13, 2, 2019)

# In BS the FX rate is the price in domestic of one unit of foreign
# In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
# DOM = USD , FOR = EUR
ccy1CCRate = 0.030  # EUR
ccy2CCRate = 0.025  # USD

strike_fx_rate = 1.250
volatility = 0.10

notional = 1000000.0

domestic_curve = DiscountCurveFlat(value_dt, ccy2CCRate)
foreign_curve = DiscountCurveFlat(value_dt, ccy1CCRate)

model = BlackScholes(volatility)

# Two examples to show that changing the notional currency and notional
# keeps the value unchanged
notional = 1000000.0

spot_fx_rates = np.arange(50, 200, 10) / 100.0


def test_european_call():
    spot_fx_rate = 1.20

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.EUROPEAN_CALL,
        notional,
        "USD",
    )
    value_european = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_european, 4) == 0.0251

    spot_fx_rate = 1.80

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.EUROPEAN_CALL,
        notional,
        "USD",
    )
    value_european = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_european, 4) == 0.5277


def test_american_call():
    spot_fx_rate = 1.20

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.AMERICAN_CALL,
        1000000,
        "USD",
    )
    value_american = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_american, 4) == 0.0255

    spot_fx_rate = 1.80

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.AMERICAN_CALL,
        1000000,
        "USD",
    )
    value_american = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_american, 4) == 0.5500


def test_european_put():
    spot_fx_rate = 1.20

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.EUROPEAN_PUT,
        notional,
        "USD",
    )
    value_european = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_european, 4) == 0.0797

    spot_fx_rate = 1.80

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.EUROPEAN_PUT,
        notional,
        "USD",
    )
    value_european = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_european, 4) == 0.0000


def test_american_put():
    spot_fx_rate = 1.20

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.AMERICAN_PUT,
        1000000,
        "USD",
    )
    value_american = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_american, 4) == 0.0798

    spot_fx_rate = 1.80

    call_option = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        "EURUSD",
        OptionTypes.AMERICAN_PUT,
        1000000,
        "USD",
    )
    value_american = call_option.value(
        value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    )["v"]

    assert round(value_american, 4) == 0.0000
