###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np

import sys
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.utils.date import Date
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXDoubleDigitalOption():

    value_date = Date(10, 4, 2020)
    expiry_date = Date(18, 9, 2020)

    forName = "EUR"
    domName = "USD"
    forCCRate = 0.03460  # EUR
    domCCRate = 0.02940  # USD

    currency_pair = forName + domName  # Always FORDOM
    spot_fx_rate = 1.20

    dom_discount_curve = DiscountCurveFlat(value_date, domCCRate)
    for_discount_curve = DiscountCurveFlat(value_date, forCCRate)

    volatility = 0.20

    notional = 1.0

    upper_strike = 1.4
    lower_strike = 1.1

    model = BlackScholes(volatility)

    double_digital_option = FXDoubleDigitalOption(
        expiry_date,
        upper_strike,
        lower_strike,
        currency_pair,
        notional,
        "USD",
    )

    spot_fx_rate = np.linspace(0.01, 2.0, 10)

    value = double_digital_option.value(
        value_date,
        spot_fx_rate,
        dom_discount_curve,
        for_discount_curve,
        model)

###############################################################################


test_FinFXDoubleDigitalOption()
testCases.compareTestCases()
