# ============================================================================
# FINANCEPY EXAMPLES - BondEmbeddedOption with Hull-White Model
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of bonds with embedded call and
# put options using the Hull-White (HW) short-rate tree.
#
# Two examples are considered:
#
#   1. A puttable bond based on a MATLAB reference example
#   2. A callable bond based on a QuantLib reference example
#
# For each example:
#
#   - a discount curve is constructed
#   - the corresponding option-free bond is valued
#   - the embedded-option bond is valued using a Hull-White tree
#   - the number of tree time steps is varied
#   - model convergence and calculation time are reported
#
# The Hull-White model is specified by:
#
#   sigma : short-rate volatility
#   a     : mean-reversion parameter
#
# BondEmbeddedOption.value() returns both the value of the bond including
# its embedded option and the corresponding option-free bond value.
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.ibor_single_curve import IborSingleCurve

from financepy.models.hw_tree import HWTree

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_embedded_option import BondEmbeddedOption
from financepy.products.rates.ibor_swap import IborSwap

set_plot_style()

# ============================================================================
# 1. PUTTABLE BOND - MATLAB REFERENCE EXAMPLE
# ============================================================================
#
# Construct an interest-rate curve from three swaps and value a bond
# containing a sequence of investor put options.
#
# The investor can put the bond back to the issuer at par on monthly dates
# beginning in January 2008. The embedded put therefore adds value relative
# to an otherwise equivalent option-free bond.
# ============================================================================

print("\n" + "=" * 90)
print("1. PUTTABLE BOND - MATLAB REFERENCE EXAMPLE")
print("=" * 90)


# ============================================================================
# 1.1 DISCOUNT CURVE
# ============================================================================
#
# Construct the discount curve from annual-pay fixed-for-floating swaps
# with maturities of one, two and three years.
# ============================================================================

settle_dt = Date(1, 1, 2007)
value_dt = settle_dt

fixed_dc_type = DayCountTypes.THIRTY_E_360
fixed_freq_type = FrequencyTypes.ANNUAL
fixed_leg_type = SwapTypes.PAY

swap1 = IborSwap(
    settle_dt,
    "1Y",
    fixed_leg_type,
    0.0350,
    fixed_freq_type,
    fixed_dc_type,
)

swap2 = IborSwap(
    settle_dt,
    "2Y",
    fixed_leg_type,
    0.0400,
    fixed_freq_type,
    fixed_dc_type,
)

swap3 = IborSwap(
    settle_dt,
    "3Y",
    fixed_leg_type,
    0.0450,
    fixed_freq_type,
    fixed_dc_type,
)

swaps = [
    swap1,
    swap2,
    swap3,
]

discount_curve = IborSingleCurve(
    value_dt,
    [],
    [],
    swaps,
)


# ============================================================================
# 1.2 OPTION-FREE BOND
# ============================================================================
#
# Define the underlying bond and value it directly from the discount curve.
# The clean price provides an independent option-free benchmark for the
# subsequent tree valuation.
# ============================================================================

issue_dt = Date(1, 1, 2004)
maturity_dt = Date(1, 1, 2010)

coupon = 0.0525
freq_type = FrequencyTypes.ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)

bond_pure_curve = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Option-free clean price from discount curve: " f"{bond_pure_curve:.6f}")


# ============================================================================
# 1.3 PUT SCHEDULE
# ============================================================================
#
# The investor can put the bond back to the issuer at 100 on monthly dates
# beginning on 1 January 2008.
# ============================================================================

call_dts = []
call_prices = []

put_dts = []
put_prices = []

put_dt = Date(1, 1, 2008)

for _ in range(24):

    put_dts.append(put_dt)
    put_prices.append(100.0)

    put_dt = put_dt.add_months(1)

call_prices = np.array(call_prices)
put_prices = np.array(put_prices)

puttable_bond = BondEmbeddedOption(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
    call_dts,
    call_prices,
    put_dts,
    put_prices,
)


# ============================================================================
# 1.4 HULL-WHITE MODEL
# ============================================================================
#
# sigma controls the volatility of the short-rate process.
#
# a controls the speed of mean reversion. A larger value causes the short
# rate to revert more rapidly toward its model-implied level.
# ============================================================================

sigma = 0.01
a = 0.10


# ============================================================================
# 1.5 TREE CONVERGENCE
# ============================================================================
#
# Revalue the puttable bond using increasingly fine Hull-White trees.
#
# The model result is unpacked into:
#
#   WITH OPTION : value including the embedded put
#   BOND PURE   : corresponding option-free value from the HW tree
#   PUT VALUE   : value added by the investor's put
#
# For a puttable bond:
#
#   put value = bond with option - option-free bond
#
# The elapsed calculation time is also reported.
# ============================================================================

print("\nHULL-WHITE TREE CONVERGENCE")
print("-" * 90)

print(f"{'STEPS':>8}" f"{'WITH OPTION':>16}" f"{'BOND PURE':>16}" f"{'PUT VALUE':>16}" f"{'TIME':>14}")

print("-" * 90)

time_steps = range(50, 1000, 50)
values = []

for num_time_steps in time_steps:

    model = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    start = time.perf_counter()

    result = puttable_bond.value(
        settle_dt,
        discount_curve,
        model,
    )

    elapsed = time.perf_counter() - start

    bond_with_option = result[0]
    bond_pure = result[1]
    put_value = bond_with_option - bond_pure

    print(
        f"{num_time_steps:8d}" f"{bond_with_option:16.6f}" f"{bond_pure:16.6f}" f"{put_value:16.6f}" f"{elapsed:14.6f}"
    )

    values.append(bond_with_option)


# Plot convergence of the puttable-bond value.
plt.figure()

plt.plot(
    list(time_steps),
    values,
    marker="o",
)

plt.title("Puttable Bond Price Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Bond Price")
plt.grid()


# ============================================================================
# 2. CALLABLE BOND - QUANTLIB REFERENCE EXAMPLE
# ============================================================================
#
# Value a callable bond using parameters based on a QuantLib reference
# example.
#
# The issuer can call the bond at par on a quarterly schedule. Because the
# call option belongs to the issuer, it reduces the value of the bond to
# the investor relative to an otherwise identical option-free bond.
# ============================================================================

print("\n" + "=" * 90)
print("2. CALLABLE BOND - QUANTLIB REFERENCE EXAMPLE")
print("=" * 90)


# ============================================================================
# 2.1 DISCOUNT CURVE
# ============================================================================
#
# Use a flat 3.5% semi-annually compounded discount curve.
# ============================================================================

value_dt = Date(16, 8, 2016)
settle_dt = value_dt.add_weekdays(3)

discount_curve = FlatDiscountCurve(
    value_dt,
    0.035,
    FrequencyTypes.SEMI_ANNUAL,
)


# ============================================================================
# 2.2 OPTION-FREE BOND
# ============================================================================

issue_dt = Date(15, 9, 2010)
maturity_dt = Date(15, 9, 2022)

coupon = 0.025
freq_type = FrequencyTypes.QUARTERLY
dc_type = DayCountTypes.ACT_ACT_ICMA

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)

bond_pure_curve = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Option-free clean price from discount curve: " f"{bond_pure_curve:.6f}")


# ============================================================================
# 2.3 CALL SCHEDULE
# ============================================================================
#
# The issuer can call the bond at par on quarterly dates beginning on
# 15 September 2016.
# ============================================================================

next_call_dt = Date(15, 9, 2016)

call_dts = [next_call_dt]
call_prices = [100.0]

for _ in range(1, 24):

    next_call_dt = next_call_dt.add_months(3)

    call_dts.append(next_call_dt)
    call_prices.append(100.0)

put_dts = []
put_prices = []

call_prices = np.array(call_prices)
put_prices = np.array(put_prices)

callable_bond = BondEmbeddedOption(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
    call_dts,
    call_prices,
    put_dts,
    put_prices,
)


# ============================================================================
# 2.4 HULL-WHITE MODEL
# ============================================================================
#
# Retain the volatility and mean-reversion assumptions from the original
# example.
#
# sigma = 0.12 is a large volatility assumption for this example, but is
# preserved here so that the results remain comparable with the original
# reference calculation.
# ============================================================================

sigma = 0.12
a = 0.03


# ============================================================================
# 2.5 TREE CONVERGENCE
# ============================================================================
#
# Compare the callable-bond value with the option-free value produced by
# the same Hull-White tree.
#
# For a callable bond:
#
#   call value = option-free bond - bond with option
#
# This difference represents the value of the issuer's embedded call.
# ============================================================================

print("\nHULL-WHITE TREE CONVERGENCE")
print("-" * 90)

print(f"{'STEPS':>8}" f"{'WITH OPTION':>16}" f"{'BOND PURE':>16}" f"{'CALL VALUE':>16}" f"{'TIME':>14}")

print("-" * 90)

time_steps = range(100, 200, 50)
values = []

for num_time_steps in time_steps:

    model = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    start = time.perf_counter()

    result = callable_bond.value(
        settle_dt,
        discount_curve,
        model,
    )

    elapsed = time.perf_counter() - start

    bond_with_option = result[0]
    bond_pure = result[1]
    call_value = bond_pure - bond_with_option

    print(
        f"{num_time_steps:8d}" f"{bond_with_option:16.6f}" f"{bond_pure:16.6f}" f"{call_value:16.6f}" f"{elapsed:14.6f}"
    )

    values.append(bond_with_option)


# Plot convergence of the callable-bond value.
plt.figure()

plt.plot(
    list(time_steps),
    values,
    marker="o",
)

plt.title("Callable Bond Price Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Bond Price")
plt.grid()

plt.show()
