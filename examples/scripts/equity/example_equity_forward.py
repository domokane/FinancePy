# ============================================================================
# FINANCEPY EXAMPLES - EquityForward
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import LongShortTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_forward import EquityForward


# ============================================================================
# 1. EQUITY FORWARD
# ============================================================================
# What this section demonstrates:
# Values an existing long equity forward and compares its contracted delivery
# price with the fair forward price implied by the current spot price,
# discount curve and dividend curve.

print("\n" + "=" * 78)
print("1. EQUITY FORWARD")
print("=" * 78)

value_dt = Date(13, 2, 2018)
expiry_dt = value_dt.add_months(12)

stock_price = 130.0
forward_price = 125.0
notional = 100.0

discount_rate = 0.05
dividend_rate = 0.02

discount_curve = FlatDiscountCurve(value_dt, discount_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

equity_forward = EquityForward(
    expiry_dt,
    forward_price,
    notional,
    LongShortTypes.LONG,
)

fair_forward = equity_forward.forward(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
)

forward_value = equity_forward.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
)

print()
print(f"{'Spot Price':>14s} {'Fair Forward':>14s} {'Forward Value':>16s}")
print("-" * 48)

print(
    f"{stock_price:14.4f}"
    f"{fair_forward:14.4f}"
    f"{forward_value:16.6f}"
)


# ============================================================================
# 2. FORWARD VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how the value of an existing long equity forward changes with the
# underlying stock price. The forward value is linear in spot under the
# deterministic discount and dividend curves used here.

print("\n" + "=" * 78)
print("2. FORWARD VALUE VERSUS STOCK PRICE")
print("=" * 78)

stock_prices = np.linspace(90.0, 160.0, 71)

forward_values = []

for stock_price_i in stock_prices:

    value = equity_forward.value(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
    )

    forward_values.append(value)

zero_value_spot = forward_price / equity_forward.forward(
    value_dt,
    1.0,
    discount_curve,
    dividend_curve,
)

fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    forward_values,
    label="Long Forward",
)

ax.axhline(
    0.0,
    linestyle="--",
    label="Zero Value",
)

ax.axvline(
    zero_value_spot,
    linestyle="--",
    label="Zero-Value Spot",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Forward Value")
ax.set_title("Long Equity Forward Value versus Stock Price")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 3. FAIR FORWARD PRICE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how the fair delivery price of a newly entered equity forward changes
# with the current stock price for fixed interest and dividend curves.

print("\n" + "=" * 78)
print("3. FAIR FORWARD PRICE VERSUS STOCK PRICE")
print("=" * 78)

fair_forward_prices = []

for stock_price_i in stock_prices:

    fair_forward_i = equity_forward.forward(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
    )

    fair_forward_prices.append(fair_forward_i)

fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    fair_forward_prices,
    label="Fair Forward Price",
)

ax.axhline(
    forward_price,
    linestyle="--",
    label="Contracted Forward Price",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Forward Price")
ax.set_title("Fair Equity Forward Price versus Stock Price")
ax.grid(True)
ax.legend()

plt.show()
