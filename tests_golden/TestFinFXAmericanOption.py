###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import numpy as np
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinFXAmericanOption():

    # There is no FXAmericanOption class. It is embedded in the FXVanillaOption
    # class. This test just compares it to the European

    valuation_date = Date(13, 2, 2018)
    expiry_date = Date(13, 2, 2019)

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

    dom_discount_curve = DiscountCurveFlat(valuation_date, ccy2CCRate)
    for_discount_curve = DiscountCurveFlat(valuation_date, ccy1CCRate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged

    testCases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spot_fx_rates = np.arange(50, 200, 10)/100.0

    for spot_fx_rate in spot_fx_rates:

        call_option = FXVanillaOption(expiry_date,
                                      strike_fx_rate,
                                      currency_pair,
                                      OptionTypes.EUROPEAN_CALL,
                                      1000000,
                                      "USD")

        valueEuropean = call_option.value(valuation_date,
                                          spot_fx_rate,
                                          dom_discount_curve,
                                          for_discount_curve,
                                          model)['v']

        call_option = FXVanillaOption(expiry_date,
                                      strike_fx_rate,
                                      "EURUSD",
                                      OptionTypes.AMERICAN_CALL,
                                      1000000,
                                      "USD")

        valueAmerican = call_option.value(valuation_date,
                                          spot_fx_rate,
                                          dom_discount_curve,
                                          for_discount_curve,
                                          model)['v']

        diff = (valueAmerican - valueEuropean)
        testCases.print(spot_fx_rate, valueEuropean, valueAmerican, diff)

    for spot_fx_rate in spot_fx_rates:

        call_option = FXVanillaOption(expiry_date,
                                      strike_fx_rate,
                                      "EURUSD",
                                      OptionTypes.EUROPEAN_PUT,
                                      1000000,
                                      "USD")

        valueEuropean = call_option.value(valuation_date,
                                          spot_fx_rate,
                                          dom_discount_curve,
                                          for_discount_curve,
                                          model)['v']

        call_option = FXVanillaOption(expiry_date,
                                      strike_fx_rate,
                                      "EURUSD",
                                      OptionTypes.AMERICAN_PUT,
                                      1000000,
                                      "USD")

        valueAmerican = call_option.value(valuation_date,
                                          spot_fx_rate,
                                          dom_discount_curve,
                                          for_discount_curve,
                                          model)['v']

        diff = (valueAmerican - valueEuropean)
        testCases.print(spot_fx_rate, valueEuropean, valueAmerican, diff)

###############################################################################


test_FinFXAmericanOption()
testCases.compareTestCases()
