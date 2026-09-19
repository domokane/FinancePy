# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.utils.global_types import TouchOptionTypes
from financepy.products.fx import FXOneTouchOption
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - FXOneTouchOption
# ============================================================================


DEBUG_FLAG = False
########################################################################################




#        print(up_type, v, v_mc)

########################################################################################




########################################################################################

# ============================================================================
# 1. FIN FX ONE TOUCH OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN FX ONE TOUCH OPTION")
print("=" * 78)

value_dt = Date(1, 1, 2016)
expiry_dt = Date(2, 7, 2016)
volatility = 0.20
barrier_level = 1.0  # H
model = BlackScholes(volatility)

domestic_rate = 0.10
foreign_rate = 0.03

num_paths = 50000
num_steps_per_year = 252 * 2

dom_curve = FlatDiscountCurve(value_dt, domestic_rate)
for_curve = FlatDiscountCurve(value_dt, foreign_rate)

spot_fx_rate = 1.050
payment_size = 1.5

print("================================= CASH ONLY")

down_types = [
    TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT,
    TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
    TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING,
]

print("TYPE", "VALUE", "VALUE_MC")

for down_type in down_types:

    option = FXOneTouchOption(expiry_dt, down_type, barrier_level, payment_size)

    v = option.value(value_dt, spot_fx_rate, dom_curve, for_curve, model)

    v_mc = option.value_mc(
        value_dt,
        spot_fx_rate,
        dom_curve,
        for_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % down_type, "%9.5f" % v, "%9.5f" % v_mc)

#        print(down_type, v, v_mc)

spot_fx_rate = 0.950
payment_size = 1.5

up_types = [
    TouchOptionTypes.UP_AND_IN_CASH_AT_HIT,
    TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
    TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING,
]

print("TYPE", "VALUE", "VALUE_MC")

for up_type in up_types:

    option = FXOneTouchOption(expiry_dt, up_type, barrier_level, payment_size)

    v = option.value(value_dt, spot_fx_rate, dom_curve, for_curve, model)

    v_mc = option.value_mc(
        value_dt,
        spot_fx_rate,
        dom_curve,
        for_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % up_type, "%9.5f" % v, "%9.5f" % v_mc)

# ============================================================================
# 2. BBG ONE TOUCH OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# Measures first-order sensitivity of value to the underlying market variable.
# Measures how delta itself changes as the underlying market variable changes.
# Measures sensitivity of value to volatility.

print("\n" + "=" * 78)
print("2. BBG ONE TOUCH OPTION")
print("=" * 78)

value_dt = Date(3, 12, 2021)
expiry_dt = Date(5, 12, 2022)
barrier_level = 1.1865  # THIS IS NUMBER OF DOLLARS PER EURO

spot_fx_rate = 1.1300  # EURUSD
volatility = 0.06075

model = BlackScholes(volatility)

for_rate = 0.00593  # EUR
dom_rate = -0.00414  # USD

num_paths = 50000
num_steps_per_year = 252

dom_curve = FlatDiscountCurve(value_dt, dom_rate)
for_curve = FlatDiscountCurve(value_dt, for_rate)

payment_size = 1000000  # EUR

option_type = TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY

option = FXOneTouchOption(expiry_dt, option_type, barrier_level, payment_size)

v = option.value(value_dt, spot_fx_rate, dom_curve, for_curve, model)

v_mc = option.value_mc(
    value_dt,
    spot_fx_rate,
    dom_curve,
    for_curve,
    model,
    num_steps_per_year,
    num_paths,
)

d = option.delta(value_dt, spot_fx_rate, dom_curve, for_curve, model)

g = option.gamma(value_dt, spot_fx_rate, dom_curve, for_curve, model)

v = option.vega(value_dt, spot_fx_rate, dom_curve, for_curve, model)

# I SHOULD GET 49.4934% OR 494,934 in EUR
# VEGA IS 68,777.26
# GAMMA IS 916,285
# DELTA IS -9560266

if DEBUG_FLAG:
    print(option_type)
    print("Value:", v)
    print("Value MC:", v_mc)
    print("Delta: ", d)
    print("Gamma:", g)
    print("Vega:", v)

