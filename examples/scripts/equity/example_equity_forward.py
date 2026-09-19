# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.

from financepy.products.equity.equity_forward import EquityForward
from financepy.utils.date import Date
from financepy.utils.global_types import LongShortTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - EquityForward
# ============================================================================


########################################################################################




########################################################################################

# ============================================================================
# 1. EQUITY FORWARD
# ============================================================================
# What this section demonstrates:
# Values the instrument using the supplied market data/model inputs. The surrounding comparison shows how the valuation responds to those assumptions.

print("\n" + "=" * 78)
print("1. EQUITY FORWARD")
print("=" * 78)

value_dt = Date(13, 2, 2018)
expiry_dt = value_dt.add_months(12)

stock_price = 130.0
forward_price = 125.0  # Locked
discount_rate = 0.05
dividend_rate = 0.02

expiry_dt = value_dt.add_months(12)
notional = 100.0

discount_curve = FlatDiscountCurve(value_dt, discount_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_rate)

equity_forward = EquityForward(expiry_dt, forward_price, notional, LongShortTypes.LONG)

print(f"{'Spot price':>14s} {'Fair forward':>14s} {'Forward value':>16s}")
print("-" * 48)

fwd_price = equity_forward.forward(value_dt, stock_price, discount_curve, dividend_curve)

fwd_value = equity_forward.value(value_dt, stock_price, discount_curve, dividend_curve)

#    print(f"{stock_price:14.4f} {fwd_price:14.4f} {fwd_value:16.6f}")
print(f"{stock_price:14.4f} {fwd_price:14.4f} {fwd_value:16.6f}")

# =============================================================================
# 2. FORWARD VALUE AS THE SPOT PRICE CHANGES
# =============================================================================
# A long equity forward becomes more valuable as the underlying spot rises.
# This sensitivity is close to linear and is clearer in a plot than from a
# single valuation point.
plot_spots = np.linspace(90.0, 160.0, 71)
plot_values = [equity_forward.value(value_dt, s, discount_curve, dividend_curve)
               for s in plot_spots]
plt.figure()
plt.plot(plot_spots, plot_values)
plt.axhline(0.0, linestyle="--")
plt.axvline(stock_price, linestyle="--", label="Example spot")
plt.xlabel("Spot price")
plt.ylabel("Forward value")
plt.title("Long equity forward value versus spot price")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
