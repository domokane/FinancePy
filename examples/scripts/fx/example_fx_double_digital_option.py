# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import numpy as np


from financepy.products.fx.fx_double_digital_option import FXDoubleDigitalOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - FXDoubleDigitalOption
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. FIN FX DOUBLE DIGITAL OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. FIN FX DOUBLE DIGITAL OPTION")
print("=" * 78)

value_dt = Date(10, 4, 2020)
expiry_dt = Date(18, 9, 2020)

for_name = "EUR"
dom_name = "USD"
for_cc_rate = 0.03460  # EUR
dom_cc_rate = 0.02940  # USD

currency_pair = for_name + dom_name  # Always FORDOM
spot_fx_rate = 1.20

domestic_curve = FlatDiscountCurve(value_dt, dom_cc_rate)
foreign_curve = FlatDiscountCurve(value_dt, for_cc_rate)

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
    value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
)

