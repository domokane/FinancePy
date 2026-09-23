# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.global_vars import ONE_MILLION
from financepy.products.rates.ois import OIS
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes

# ============================================================================
# FINANCEPY EXAMPLES - Ois
# ============================================================================


########################################################################################


########################################################################################

# ============================================================================
# 1. FIN FIXED OIS
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. FIN FIXED OIS")
print("=" * 78)

effective_dt = Date(30, 11, 2018)
end_dt = Date(30, 11, 2023)

end_dt = effective_dt.add_months(60)
ois_rate = 0.04
fixed_leg_type = SwapTypes.PAY
fixed_freq_type = FrequencyTypes.ANNUAL
fixed_day_count = DayCountTypes.ACT_360
float_freq_type = FrequencyTypes.ANNUAL
float_day_count = DayCountTypes.ACT_360
float_spread = 0.0
notional = ONE_MILLION
payment_lag = 1

ois = OIS(
    effective_dt,
    end_dt,
    fixed_leg_type,
    ois_rate,
    fixed_freq_type,
    fixed_day_count,
    notional,
    payment_lag,
    float_spread,
    float_freq_type,
    float_day_count,
)

#    print(ois)

value_dt = effective_dt
market_rate = 0.05
ois_curve = FlatDiscountCurve(value_dt, market_rate, FrequencyTypes.ANNUAL)

v = ois.value(effective_dt, ois_curve)

#    print(v)

#    ois._fixed_leg.print_valuation()
#    ois._float_leg.print_valuation()

print("LABEL", "VALUE")
print("SWAP_VALUE", v)
