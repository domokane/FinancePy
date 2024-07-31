###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import sys
sys.path.append("..")

import numpy as np

from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode


test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXDoubleDigitalOption():

    value_dt = Date(10, 4, 2020)
    expiry_dt = Date(18, 9, 2020)

    for_name = "EUR"
    dom_name = "USD"
    for_cc_rate = 0.03460  # EUR
    dom_cc_rate = 0.02940  # USD

    currency_pair = for_name + dom_name  # Always FORDOM
    spot_fx_rate = 1.20

    domestic_curve = DiscountCurveFlat(value_dt, dom_cc_rate)
    foreign_curve = DiscountCurveFlat(value_dt, for_cc_rate)

    volatility = 0.20

    notional = 1.0

    upper_strike = 1.4
    lower_strike = 1.1

    model = BlackScholes(volatility)

    double_digital_option = FXDoubleDigitalOption(
        expiry_dt,
        upper_strike,
        lower_strike,
        currency_pair,
        notional,
        "USD",
    )

    spot_fx_rate = np.linspace(0.01, 2.0, 10)

    value = double_digital_option.value(
        value_dt,
        spot_fx_rate,
        domestic_curve,
        foreign_curve,
        model)

###############################################################################


test_FinFXDoubleDigitalOption()
test_cases.compareTestCases()
