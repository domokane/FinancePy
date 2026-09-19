# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.



from financepy.products.equity.equity_one_touch_option import (
    EquityOneTouchOption,
)
from financepy.utils.global_types import TouchOptionTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date

# ============================================================================
# FINANCEPY EXAMPLES - EquityOneTouchOption
# ============================================================================

########################################################################################




########################################################################################

# ============================================================================
# 1. EQUITY ONE TOUCH OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. EQUITY ONE TOUCH OPTION")
print("=" * 78)

value_dt = Date(1, 1, 2016)
expiry_dt = Date(2, 7, 2016)
interest_rate = 0.10
volatility = 0.20
barrier_level = 100.0  # H
model = BlackScholes(volatility)
dividend_yield = 0.03
num_paths = 10000
num_steps_per_year = 252

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

stock_price = 105.0
payment_size = 15.0

print("================================= CASH ONLY")

down_types = [
    TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT,
    TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
    TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING,
]

print("TYPE", "VALUE", "VALUE_MC")

for down_type in down_types:

    option = EquityOneTouchOption(expiry_dt, down_type, barrier_level, payment_size)

    v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % down_type, "%9.5f" % v, "%9.5f" % v_mc)

stock_price = 95.0
payment_size = 15.0

up_types = [
    TouchOptionTypes.UP_AND_IN_CASH_AT_HIT,
    TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
    TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING,
]

print("TYPE", "VALUE", "VALUE_MC")

for up_type in up_types:

    option = EquityOneTouchOption(expiry_dt, up_type, barrier_level, payment_size)

    v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % up_type, "%9.5f" % v, "%9.5f" % v_mc)

stock_price = 105.0

print("================= ASSET ONLY")

down_types = [
    TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT,
    TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY,
    TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING,
]

print("TYPE", "VALUE", "VALUE_MC")
for down_type in down_types:

    option = EquityOneTouchOption(expiry_dt, down_type, barrier_level)

    v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % down_type, "%9.5f" % v, "%9.5f" % v_mc)

stock_price = 95.0

up_types = [
    TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT,
    TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY,
    TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING,
]

for up_type in up_types:

    option = EquityOneTouchOption(expiry_dt, up_type, barrier_level)

    v = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

    v_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_steps_per_year,
        num_paths,
    )

    print("%60s " % up_type, "%9.5f" % v, "%9.5f" % v_mc)

