# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

import numpy as np

from financepy.products.fx import FXDoubleOneTouchOption
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date
from financepy.utils.global_types import DoubleBarrierTypes

# ============================================================================
# FINANCEPY EXAMPLES - FXDoubleOneTouchOption
# ============================================================================


DEBUG_FLAG = False
########################################################################################




##################################################################################




########################################################################################

# ============================================================================
# 1. DOUBLE ONE TOUCH SYMMETRIC BARRIERS
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The original example contains an accuracy/consistency assertion; the surrounding values show what is being checked numerically.

print("\n" + "=" * 78)
print("1. DOUBLE ONE TOUCH SYMMETRIC BARRIERS")
print("=" * 78)

"""Regression: a zero Fourier term must not terminate the DNT series."""

value_dt = Date(
    1,
    1,
    2026,
)
expiry_dt = value_dt.add_days(
    30,
)

sigma = 0.20

domestic_curve = FlatDiscountCurve(
    value_dt,
    0.5 * sigma * sigma,
)

foreign_curve = FlatDiscountCurve(
    value_dt,
    0.0,
)

model = BlackScholes(
    sigma,
)

spot = 100.0
lower = 80.0
upper = 125.0
payment = 1.0

knock_out = FXDoubleOneTouchOption(
    expiry_dt,
    DoubleBarrierTypes.KNOCK_OUT,
    lower,
    upper,
    payment,
)

knock_in = FXDoubleOneTouchOption(
    expiry_dt,
    DoubleBarrierTypes.KNOCK_IN,
    lower,
    upper,
    payment,
)

value_out = knock_out.value(
    value_dt,
    spot,
    domestic_curve,
    foreign_curve,
    model,
)

value_in = knock_in.value(
    value_dt,
    spot,
    domestic_curve,
    foreign_curve,
    model,
)

discounted_payment = payment * domestic_curve.df(expiry_dt)

# Both legs are discounted indicator payoffs and must therefore
# lie between zero and the discounted payment.
assert 0.0 <= value_out <= discounted_payment
assert 0.0 <= value_in <= discounted_payment

# A touch and no-touch partition all possible paths.
assert np.isclose(
    value_in + value_out,
    discounted_payment,
    rtol=1e-12,
    atol=1e-12,
)

expected_dnt = 0.9981587590142644

assert np.isclose(
    value_out,
    expected_dnt,
    rtol=1e-12,
    atol=1e-12,
)

assert 0.0 <= value_out <= discounted_payment
assert 0.0 <= value_in <= discounted_payment

assert np.isclose(
    value_in + value_out,
    discounted_payment,
    rtol=1e-12,
    atol=1e-12,
)

# ============================================================================
# 2. FIN FX DOUBLE BARRIER ONE TOUCH OPTION
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("2. FIN FX DOUBLE BARRIER ONE TOUCH OPTION")
print("=" * 78)

value_dt = Date(1, 1, 2016)
expiry_dt = value_dt.add_months(3)

spot_fx_rate = 100.0

domestic_rate = 0.05
foreign_rate = 0.05 - 0.03

num_paths = 10000
num_steps_per_year = 252 * 2

dom_curve = FlatDiscountCurve(value_dt, domestic_rate)
for_curve = FlatDiscountCurve(value_dt, foreign_rate)

payment_size = 10.0

lower_barrier_fx_rate = [80.0, 85.0, 90.0, 95.0]
upper_barrier_fx_rate = [120.0, 115.0, 110.0, 105.0]

for option_type in [DoubleBarrierTypes.KNOCK_OUT, DoubleBarrierTypes.KNOCK_IN]:

    for i in range(0, 4):

        k1 = lower_barrier_fx_rate[i]
        k2 = upper_barrier_fx_rate[i]

        option = FXDoubleOneTouchOption(expiry_dt, option_type, k1, k2, payment_size)

        for sigma in [0.1, 0.2, 0.3, 0.5]:

            volatility = sigma
            model = BlackScholes(volatility)

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

            #                print(option_type, k1, k2, sigma, v, v_mc)

            print("=================================")
            print("OPT_TYPE", "L", "U", "SIGMA", "VALUE", "VALUE_MC")
            print(option_type, k1, k2, sigma, v, v_mc)

