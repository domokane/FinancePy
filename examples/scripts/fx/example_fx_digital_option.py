# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.utils.global_types import OptionTypes
from financepy.products.fx.fx_digital_option import FXDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - FXDigitalOption
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. FIN FX DIGITAL OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. FIN FX DIGITAL OPTION")
print("=" * 78)

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

notional = 1.0

domestic_curve = FlatDiscountCurve(value_dt, ccy2_cc_rate)
foreign_curve = FlatDiscountCurve(value_dt, ccy1_cc_rate)

model = BlackScholes(volatility)

digital_option = FXDigitalOption(
    expiry_dt,
    strike_fx_rate,
    currency_pair,
    OptionTypes.DIGITAL_CALL,
    notional,
    "USD",
)

spot_fx_rate = np.linspace(0.01, 2.0, 10)

value = digital_option.value(value_dt, spot_fx_rate, domestic_curve, foreign_curve, model)

