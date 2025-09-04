# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date

from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_fx_american_option():

    # There is no FXAmericanOption class. It is embedded in the FXVanillaOption
    # class. This test just compares it to the European

    value_dt = Date(13, 2, 2018)
    expiry_dt = Date(13, 2, 2019)

    # In BS the FX rate is the price in domestic of one unit of foreign
    # In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
    # DOM = USD , FOR = EUR
    ccy1 = "EUR"
    ccy2 = "USD"
    ccy1_cc_rate = 0.030  # EUR
    ccy2_cc_rate = 0.025  # USD

    currency_pair = ccy1 + ccy2  # Always ccy1ccy2
    spot_fx_rate = 1.20
    strike_fx_rate = 1.250
    volatility = 0.10

    domestic_curve = DiscountCurveFlat(value_dt, ccy2_cc_rate)
    foreign_curve = DiscountCurveFlat(value_dt, ccy1_cc_rate)

    model = BlackScholes(volatility)

    # Two examples to show that changing the notional currency and notional
    # keeps the value unchanged

    test_cases.header("SPOT FX RATE", "VALUE_BS", "VOL_IN", "IMPLD_VOL")

    spot_fx_rates = np.arange(50, 200, 20) / 100.0

    for spot_fx_rate in spot_fx_rates:

        call_option = FXVanillaOption(
            expiry_dt,
            strike_fx_rate,
            currency_pair,
            OptionTypes.EUROPEAN_CALL,
            1000000,
            "USD",
        )

        value_european = call_option.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )["v"]

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

        diff = value_american - value_european
        test_cases.print(spot_fx_rate, value_european, value_american, diff)

    for spot_fx_rate in spot_fx_rates:

        call_option = FXVanillaOption(
            expiry_dt,
            strike_fx_rate,
            "EURUSD",
            OptionTypes.EUROPEAN_PUT,
            1000000,
            "USD",
        )

        value_european = call_option.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )["v"]

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

        diff = value_american - value_european
        test_cases.print(spot_fx_rate, value_european, value_american, diff)


########################################################################################

test_fin_fx_american_option()
test_cases.compare_test_cases()
