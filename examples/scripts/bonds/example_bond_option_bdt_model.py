# ============================================================================
# FINANCEPY EXAMPLES - BondOption BDT Model
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of European and American options
# on fixed-rate bonds using the Black-Derman-Toy (BDT) interest-rate tree.
#
# Three examples are considered:
#
#   1. Near-zero-volatility convergence against intrinsic values
#   2. European and American call/put values for different strikes
#   3. Tree convergence of European and American option values
#
# The examples illustrate the effects of:
#
#   - option strike
#   - exercise style
#   - interest-rate volatility
#   - number of tree time steps
#
# Set PLOT_GRAPHS to True to display the convergence plots.
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

from financepy.models.bdt_tree import BDTTree

from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_option import BondOption

set_plot_style()

# ============================================================================
# 1. NEAR-ZERO-VOLATILITY CONVERGENCE
# ============================================================================
#
# At very low interest-rate volatility, the option values should approach
# their deterministic values.
#
# European option values can therefore be compared with the discounted
# intrinsic value calculated from the forward clean bond price.
#
# American options are also calculated to show whether early exercise adds
# value under these assumptions.
# ============================================================================

print("\n" + "=" * 100)
print("1. BOND OPTION - NEAR-ZERO-VOLATILITY CONVERGENCE")
print("=" * 100)


# ============================================================================
# 1.1 BOND AND DISCOUNT CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

rate = 0.05

discount_curve = FlatDiscountCurve(
    settle_dt,
    rate,
    FrequencyTypes.ANNUAL,
)

issue_dt = Date(1, 9, 2015)
maturity_dt = Date(1, 9, 2025)

coupon = 0.06
freq_type = FrequencyTypes.ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)


# ============================================================================
# 1.2 OPTION EXPIRY AND REFERENCE VALUES
# ============================================================================

expiry_dt = settle_dt.add_tenor("18m")

df_expiry = discount_curve.df(expiry_dt)

spot_clean_value = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

fwd_clean_value = bond.clean_price_from_discount_curve(
    expiry_dt,
    discount_curve,
)

print(f"Spot clean bond price    : {spot_clean_value:.6f}")
print(f"Forward clean bond price : {fwd_clean_value:.6f}")
print(f"Expiry discount factor   : {df_expiry:.6f}")


# ============================================================================
# 1.3 CONVERGENCE TEST
# ============================================================================
#
# A very small volatility is used so that the BDT tree is effectively
# approaching a deterministic interest-rate model.
# ============================================================================

strike_prices = [
    90.0,
    100.0,
    110.0,
    120.0,
]

num_time_steps = range(100, 1000, 200)

sigma = 1.0e-7

print("\n" + "-" * 118)

print(
    f"{'STRIKE':>8}"
    f"{'STEPS':>8}"
    f"{'CALL INT':>13}"
    f"{'CALL INT PV':>14}"
    f"{'CALL EUR':>13}"
    f"{'CALL AMER':>13}"
    f"{'PUT INT':>13}"
    f"{'PUT INT PV':>14}"
    f"{'PUT EUR':>13}"
    f"{'PUT AMER':>13}"
)

print("-" * 118)

for strike_price in strike_prices:

    call_intrinsic = max(
        spot_clean_value - strike_price,
        0.0,
    )

    put_intrinsic = max(
        strike_price - spot_clean_value,
        0.0,
    )

    call_intrinsic_pv = (
        max(
            fwd_clean_value - strike_price,
            0.0,
        )
        * df_expiry
    )

    put_intrinsic_pv = (
        max(
            strike_price - fwd_clean_value,
            0.0,
        )
        * df_expiry
    )

    for num_steps in num_time_steps:

        model = BDTTree(
            sigma,
            num_steps,
        )

        euro_call = BondOption(
            bond,
            expiry_dt,
            strike_price,
            OptionTypes.EUROPEAN_CALL,
        )

        amer_call = BondOption(
            bond,
            expiry_dt,
            strike_price,
            OptionTypes.AMERICAN_CALL,
        )

        euro_put = BondOption(
            bond,
            expiry_dt,
            strike_price,
            OptionTypes.EUROPEAN_PUT,
        )

        amer_put = BondOption(
            bond,
            expiry_dt,
            strike_price,
            OptionTypes.AMERICAN_PUT,
        )

        value_euro_call = euro_call.value(
            settle_dt,
            discount_curve,
            model,
        )

        value_amer_call = amer_call.value(
            settle_dt,
            discount_curve,
            model,
        )

        value_euro_put = euro_put.value(
            settle_dt,
            discount_curve,
            model,
        )

        value_amer_put = amer_put.value(
            settle_dt,
            discount_curve,
            model,
        )

        print(
            f"{strike_price:8.2f}"
            f"{num_steps:8d}"
            f"{call_intrinsic:13.6f}"
            f"{call_intrinsic_pv:14.6f}"
            f"{value_euro_call:13.6f}"
            f"{value_amer_call:13.6f}"
            f"{put_intrinsic:13.6f}"
            f"{put_intrinsic_pv:14.6f}"
            f"{value_euro_put:13.6f}"
            f"{value_amer_put:13.6f}"
        )


# ============================================================================
# 2. BOND OPTION VALUES BY STRIKE
# ============================================================================
#
# Value European and American calls and puts for several strikes using a
# BDT tree with a fixed number of time steps.
#
# The discount curve is constructed directly from discount factors.
# ============================================================================

print("\n" + "=" * 100)
print("2. BOND OPTION VALUES BY STRIKE")
print("=" * 100)


# ============================================================================
# 2.1 BOND
# ============================================================================

settle_dt = Date(1, 12, 2019)
issue_dt = Date(1, 12, 2018)
maturity_dt = settle_dt.add_tenor("10Y")

coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)


# ============================================================================
# 2.2 DISCOUNT CURVE
# ============================================================================
#
# Construct a discount curve corresponding to a continuously compounded
# flat 5% rate.
# ============================================================================

t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR

times = np.linspace(
    0.0,
    t_mat,
    20,
)

dates = settle_dt.add_years(times)

dfs = np.exp(-0.05 * times)

discount_curve = DiscountCurve(
    settle_dt,
    dates,
    dfs,
)


# ============================================================================
# 2.3 OPTION INPUTS
# ============================================================================

expiry_dt = settle_dt.add_tenor("18m")

strikes = [
    80.0,
    100.0,
    120.0,
]

sigma = 0.20
num_time_steps = 100

bond_price = bond.dirty_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Dirty bond price : {bond_price:.6f}")
print(f"BDT volatility   : {sigma:.6f}")
print(f"BDT time steps   : {num_time_steps}")


# ============================================================================
# 2.4 OPTION VALUES
# ============================================================================
#
# Use the same BDT model assumptions for each option so that the effects of
# strike and exercise style can be compared directly.
# ============================================================================

print("\n" + "-" * 76)

print(f"{'OPTION TYPE':<24}" f"{'STRIKE':>12}" f"{'VALUE':>16}")

print("-" * 76)

option_types = [
    ("EUROPEAN CALL", OptionTypes.EUROPEAN_CALL),
    ("AMERICAN CALL", OptionTypes.AMERICAN_CALL),
    ("EUROPEAN PUT", OptionTypes.EUROPEAN_PUT),
    ("AMERICAN PUT", OptionTypes.AMERICAN_PUT),
]

for option_label, option_type in option_types:

    for strike_price in strikes:

        bond_option = BondOption(
            bond,
            expiry_dt,
            strike_price,
            option_type,
        )

        model = BDTTree(
            sigma,
            num_time_steps,
        )

        option_value = bond_option.value(
            settle_dt,
            discount_curve,
            model,
        )

        print(f"{option_label:<24}" f"{strike_price:12.2f}" f"{option_value:16.6f}")


# ============================================================================
# 3. AMERICAN AND EUROPEAN OPTION CONVERGENCE
# ============================================================================
#
# Examine convergence of European and American call and put values as the
# number of BDT tree time steps is increased.
#
# All four options have the same strike and expiry, making the effect of
# exercise style directly observable.
# ============================================================================

print("\n" + "=" * 100)
print("3. AMERICAN AND EUROPEAN OPTION CONVERGENCE")
print("=" * 100)


# ============================================================================
# 3.1 BOND AND CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

discount_curve = FlatDiscountCurve(
    settle_dt,
    0.05,
    FrequencyTypes.CONTINUOUS,
)

issue_dt = Date(1, 9, 2014)
maturity_dt = Date(1, 9, 2025)

coupon = 0.05
freq_type = FrequencyTypes.ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    dc_type,
)

expiry_dt = settle_dt.add_tenor("18m")

spot_value = bond.dirty_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Dirty bond price : {spot_value:.6f}")


# ============================================================================
# 3.2 OPTION AND MODEL INPUTS
# ============================================================================

strike_price = 100.0
sigma = 0.20

print(f"Strike           : {strike_price:.2f}")
print(f"BDT volatility   : {sigma:.6f}")


# ============================================================================
# 3.3 CONVERGENCE
# ============================================================================
#
# The original example used:
#
#     range(100, 100, 1)
#
# which is empty. Use a genuine convergence range here.
# ============================================================================

num_steps_vector = range(
    100,
    401,
    25,
)

values_euro_call = []
values_amer_call = []
values_euro_put = []
values_amer_put = []

print("\n" + "-" * 92)

print(f"{'STEPS':>8}" f"{'EUR CALL':>16}" f"{'AMER CALL':>16}" f"{'EUR PUT':>16}" f"{'AMER PUT':>16}" f"{'TIME':>16}")

print("-" * 92)

for num_steps in num_steps_vector:

    model = BDTTree(
        sigma,
        num_steps,
    )

    start = time.perf_counter()

    euro_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    amer_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_CALL,
    )

    euro_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_PUT,
    )

    amer_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_PUT,
    )

    value_euro_call = euro_call.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_amer_call = amer_call.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_euro_put = euro_put.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_amer_put = amer_put.value(
        settle_dt,
        discount_curve,
        model,
    )

    elapsed = time.perf_counter() - start

    print(
        f"{num_steps:8d}"
        f"{value_euro_call:16.6f}"
        f"{value_amer_call:16.6f}"
        f"{value_euro_put:16.6f}"
        f"{value_amer_put:16.6f}"
        f"{elapsed:16.6f}"
    )

    values_euro_call.append(value_euro_call)
    values_amer_call.append(value_amer_call)
    values_euro_put.append(value_euro_put)
    values_amer_put.append(value_amer_put)


# ============================================================================
# 3.4 CONVERGENCE PLOTS
# ============================================================================
#
# Plot the convergence of each option value as the number of BDT tree
# time steps is increased.
# ============================================================================

steps = list(num_steps_vector)

plt.figure()
plt.plot(
    steps,
    values_euro_call,
    label="European Call",
)
plt.title("European Call Convergence - BDT Model")
plt.xlabel("Number of BDT Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_call,
    label="American Call",
)
plt.title("American Call Convergence - BDT Model")
plt.xlabel("Number of BDT Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_euro_put,
    label="European Put",
)
plt.title("European Put Convergence - BDT Model")
plt.xlabel("Number of BDT Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_put,
    label="American Put",
)
plt.title("American Put Convergence - BDT Model")
plt.xlabel("Number of BDT Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.show()
