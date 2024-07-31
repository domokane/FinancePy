###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import sys
sys.path.append("..")

import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.utils.date import Date
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXDigitalOption():

    # Not exactly T=1.0 but close so don't exact exact agreement
    # (in fact I do not get exact agreement even if I do set T=1.0)
    value_dt = Date(13, 2, 2018)
    expiry_dt = Date(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1CCRate = 0.030  # EUR
    ccy2CCRate = 0.025  # USD

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 1.20
    strike_fx_rate = 1.250
    volatility = 0.10

    notional = 1.0

    domestic_curve = DiscountCurveFlat(value_dt, ccy2CCRate)
    foreign_curve = DiscountCurveFlat(value_dt, ccy1CCRate)

    model = BlackScholes(volatility)

    digital_option = FXDigitalOption(expiry_dt,
                                     strike_fx_rate,
                                     currency_pair,
                                     OptionTypes.DIGITAL_CALL,
                                     notional,
                                     "USD")

    spot_fx_rate = np.linspace(0.01, 2.0, 10)

    value = digital_option.value(value_dt,
                                 spot_fx_rate,
                                 domestic_curve,
                                 foreign_curve,
                                 model)

###############################################################################


test_FinFXDigitalOption()
test_cases.compareTestCases()
