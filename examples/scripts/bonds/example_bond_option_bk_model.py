# ============================================================================
# FINANCEPY EXAMPLES - BondOption with Black-Karasinski Model
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of European and American options
# on fixed-rate bonds using the Black-Karasinski (BK) short-rate tree.
#
# Four examples are considered:
#
#   1. Near-zero-volatility convergence
#   2. Bond option values for different strikes and model parameters
#   3. European and American option convergence
#   4. Detailed tree convergence with plots
#
# The examples illustrate the effects of:
#
#   - option strike
#   - exercise style
#   - short-rate volatility
#   - mean reversion
#   - number of tree time steps
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.bk_tree import BKTree
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# 1. BOND OPTION - NEAR-ZERO-VOLATILITY CONVERGENCE
# ============================================================================
#
# At very low short-rate volatility, the option values should approach
# their deterministic values.
#
# European option values can therefore be compared with discounted
# intrinsic values calculated from the forward clean bond price.
# ============================================================================

print("\n" + "=" * 100)
print("1. BOND OPTION - NEAR-ZERO-VOLATILITY CONVERGENCE")
print("=" * 100)


# ============================================================================
# 1.1 BOND AND DISCOUNT CURVE
# ============================================================================

settle_dt = Date(1, 9, 2019)

rate = 0.05

discount_curve = FlatDiscountCurve(
    settle_dt,
    rate,
    FrequencyTypes.ANNUAL,
)

issue_dt = Date(1, 9, 2014)
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

expiry_dt = Date(1, 12, 2021)

df_expiry = discount_curve.df(expiry_dt)

fwd_clean_value = bond.clean_price_from_discount_curve(
    expiry_dt,
    discount_curve,
)

fwd_full_value = bond.dirty_price_from_discount_curve(
    expiry_dt,
    discount_curve,
)

spot_clean_value = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Spot clean bond price    : {spot_clean_value:.6f}")
print(f"Forward clean bond price : {fwd_clean_value:.6f}")
print(f"Forward dirty bond price : {fwd_full_value:.6f}")
print(f"Expiry discount factor   : {df_expiry:.6f}")


# ============================================================================
# 1.3 CONVERGENCE TEST
# ============================================================================

num_time_steps = range(100, 200, 100)

strike_prices = [
    90.0,
    100.0,
    110.0,
]

sigma = 1.0e-7
a = 0.10

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

        model = BKTree(
            sigma,
            a,
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
# Value European and American calls and puts for several strikes using the
# Black-Karasinski model.
#
# Two sets of model parameters are used to demonstrate the effect of
# volatility and mean reversion on option values.
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
    90.0,
    100.0,
    110.0,
    120.0,
]

bond_price = bond.dirty_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Dirty bond price : {bond_price:.6f}")


# ============================================================================
# 2.4 OPTION VALUES
# ============================================================================
#
# Each row reports:
#
#   - exercise style
#   - strike
#   - BK volatility
#   - BK mean reversion
#   - option value
# ============================================================================

print("\n" + "-" * 88)

print(f"{'OPTION TYPE':<22}" f"{'STRIKE':>10}" f"{'SIGMA':>12}" f"{'MEAN REV':>12}" f"{'VALUE':>16}")

print("-" * 88)

option_types = [
    ("EUROPEAN CALL", OptionTypes.EUROPEAN_CALL),
    ("AMERICAN CALL", OptionTypes.AMERICAN_CALL),
    ("EUROPEAN PUT", OptionTypes.EUROPEAN_PUT),
    ("AMERICAN PUT", OptionTypes.AMERICAN_PUT),
]

parameter_sets = [
    (0.01, 0.10),
    (0.20, 0.05),
]

num_time_steps = 20

for option_label, option_type in option_types:

    for sigma, a in parameter_sets:

        for strike_price in strikes:

            bond_option = BondOption(
                bond,
                expiry_dt,
                strike_price,
                option_type,
            )

            model = BKTree(
                sigma,
                a,
                num_time_steps,
            )

            option_value = bond_option.value(
                settle_dt,
                discount_curve,
                model,
            )

            print(f"{option_label:<22}" f"{strike_price:10.2f}" f"{sigma:12.6f}" f"{a:12.6f}" f"{option_value:16.6f}")


# ============================================================================
# 3. BOND OPTION - AMERICAN CONVERGENCE
# ============================================================================
#
# Compare European and American put and call values as the number of
# Black-Karasinski tree time steps is increased.
# ============================================================================

print("\n" + "=" * 100)
print("3. BOND OPTION - AMERICAN CONVERGENCE")
print("=" * 100)


# ============================================================================
# 3.1 BOND AND CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

discount_curve = FlatDiscountCurve(
    settle_dt,
    0.05,
)

issue_dt = Date(1, 9, 2016)
maturity_dt = Date(1, 9, 2025)

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

expiry_dt = Date(1, 12, 2020)
strike_price = 100.0

sigma = 0.20
a = 0.10


# ============================================================================
# 3.2 CONVERGENCE
# ============================================================================

time_steps = range(
    20,
    100,
    20,
)

print("\n" + "-" * 94)

print(f"{'STEPS':>8}" f"{'AMER PUT':>16}" f"{'EUR PUT':>16}" f"{'AMER CALL':>16}" f"{'EUR CALL':>16}" f"{'TIME':>16}")

print("-" * 94)

for num_time_steps in time_steps:

    model = BKTree(
        sigma,
        a,
        num_time_steps,
    )

    start = time.perf_counter()

    amer_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_PUT,
    )

    euro_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_PUT,
    )

    amer_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_CALL,
    )

    euro_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    value_amer_put = amer_put.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_euro_put = euro_put.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_amer_call = amer_call.value(
        settle_dt,
        discount_curve,
        model,
    )

    value_euro_call = euro_call.value(
        settle_dt,
        discount_curve,
        model,
    )

    elapsed = time.perf_counter() - start

    print(
        f"{num_time_steps:8d}"
        f"{value_amer_put:16.6f}"
        f"{value_euro_put:16.6f}"
        f"{value_amer_call:16.6f}"
        f"{value_euro_call:16.6f}"
        f"{elapsed:16.6f}"
    )


# ============================================================================
# 4. BOND OPTION - DETAILED TREE CONVERGENCE
# ============================================================================
#
# Examine convergence of European and American call and put values using
# increasingly fine Black-Karasinski trees.
#
# The resulting convergence series are plotted automatically.
# ============================================================================

print("\n" + "=" * 100)
print("4. BOND OPTION - DETAILED TREE CONVERGENCE")
print("=" * 100)


# ============================================================================
# 4.1 BOND AND CURVE
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
# 4.2 OPTION AND MODEL INPUTS
# ============================================================================

strike_price = 100.0

sigma = 0.20
a = 0.10

print(f"Strike           : {strike_price:.2f}")
print(f"BK volatility    : {sigma:.6f}")
print(f"BK mean reversion: {a:.6f}")


# ============================================================================
# 4.3 REFERENCE OPTION VALUE
# ============================================================================

reference_model = BKTree(
    sigma,
    a,
    100,
)

reference_option = BondOption(
    bond,
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
)

reference_value = reference_option.value(
    settle_dt,
    discount_curve,
    reference_model,
)

print(f"European call value " f"(100 steps): {reference_value:.6f}")


# ============================================================================
# 4.4 TREE CONVERGENCE
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

    model = BKTree(
        sigma,
        a,
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
# 4.5 CONVERGENCE PLOTS
# ============================================================================
#
# Plot all convergence results. These plots are always generated when the
# example is run.
# ============================================================================

steps = list(num_steps_vector)

plt.figure()
plt.plot(
    steps,
    values_euro_call,
    label="European Call",
)
plt.title("European Call Convergence - BK Model")
plt.xlabel("Number of BK Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_call,
    label="American Call",
)
plt.title("American Call Convergence - BK Model")
plt.xlabel("Number of BK Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_euro_put,
    label="European Put",
)
plt.title("European Put Convergence - BK Model")
plt.xlabel("Number of BK Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_put,
    label="American Put",
)
plt.title("American Put Convergence - BK Model")
plt.xlabel("Number of BK Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.show()
