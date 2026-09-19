# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - FXVanillaOption
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. FIN FX OPTION SABR
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN FX OPTION SABR")
print("=" * 78)

value_dt = Date(13, 2, 2018)
expiry_dt = Date(13, 2, 2019)

# In BS the FX rate is the price in domestic of one unit of foreign
# In case of EURUSD = 1.3 the domestic currency is USD and foreign is EUR
# DOM = USD , FOR = EUR
ccy1_cc_rate = 0.030  # EUR
ccy2_cc_rate = 0.025  # USD

spot_fx_rate = 1.20
strike_fx_rate = 1.250
volatility = 0.10

notional = 1000000.0

domestic_curve = FlatDiscountCurve(value_dt, ccy2_cc_rate)
foreign_curve = FlatDiscountCurve(value_dt, ccy1_cc_rate)

model = BlackScholes(volatility)

# Two examples to show that changing the notional currency and notional
# keeps the value unchanged
notional = 1000000.0

spot_fx_rates = np.arange(50, 200, 10) / 100.0

print("OPTION", "FX_RATE", "VALUE_BS", "VOL_IN", "DIFF")

for spot_fx_rate in spot_fx_rates:

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

    print(
        "CALL:",
        "%9.6f" % spot_fx_rate,
        "%9.7f" % value_european,
        "%9.7f" % value_american,
        "%9.7f" % diff,
    )

print("OPTION", "FX_RATE", "VALUE_BS", "VOL_IN", "DIFF")

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
    print(
        "PUT:",
        "%9.6f" % spot_fx_rate,
        "%9.7f" % value_european,
        "%9.7f" % value_american,
        "%9.7f" % diff,
    )

