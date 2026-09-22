# ============================================================================
# FINANCEPY EXAMPLES - BondOption with Hull-White Model
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of European and American options
# on fixed-rate bonds using the Hull-White (HW) short-rate model.
#
# The examples cover:
#
#   1. Comparison of Jamshidian and tree-based European option valuation
#   2. Near-zero-volatility convergence
#   3. Bond option values across different strikes
#   4. European option calculation convergence
#   5. American versus European option convergence
#   6. Detailed tree convergence with plots
#
# The examples illustrate the effects of:
#
#   - option strike
#   - exercise style
#   - short-rate volatility
#   - mean reversion
#   - number of tree time steps
#   - European option calculation method
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.hw_tree import HWTree, HWEuropeanCalcTypes
from financepy.products.bonds.bond import Bond
from financepy.products.bonds.bond_option import BondOption
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import ExerciseTypes
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_vars import G_DAYS_IN_YEAR
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# 1. BOND OPTION - JAMSHIDIAN VERSUS TREE
# ============================================================================
#
# Compare two approaches for valuing a European option on a coupon bond
# under the Hull-White model:
#
#   - Jamshidian decomposition
#   - Hull-White interest-rate tree
#
# Both methods return call and put option values. These are unpacked and
# compared directly.
# ============================================================================

print("\n" + "=" * 100)
print("1. BOND OPTION - JAMSHIDIAN VERSUS TREE")
print("=" * 100)


# ============================================================================
# 1.1 BOND AND DISCOUNT CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

rate = 0.05
fixed_dc_type = DayCountTypes.THIRTY_360_BOND
fixed_freq_type = FrequencyTypes.SEMI_ANNUAL

discount_curve = FlatDiscountCurve(
    settle_dt,
    rate,
    fixed_freq_type,
    fixed_dc_type,
)

issue_dt = Date(1, 12, 2018)
expiry_dt = settle_dt.add_tenor("18m")
maturity_dt = settle_dt.add_tenor("10Y")

coupon = 0.05
freq_type = FrequencyTypes.SEMI_ANNUAL
accrual_type = DayCountTypes.THIRTY_360_BOND

bond = Bond(
    issue_dt,
    maturity_dt,
    coupon,
    freq_type,
    accrual_type,
)

strike_price = 100.0
face = 100.0


# ============================================================================
# 1.2 COUPON CASH FLOWS
# ============================================================================

coupon_times = []
coupon_flows = []

cpn = bond.cpn / bond.freq

num_flows = len(bond.cpn_dts)

for i in range(num_flows):

    pcd = bond.cpn_dts[i - 1]
    ncd = bond.cpn_dts[i]

    if ncd > settle_dt:

        if len(coupon_times) == 0:

            flow_time = (pcd - settle_dt) / G_DAYS_IN_YEAR

            coupon_times.append(flow_time)
            coupon_flows.append(cpn)

        flow_time = (ncd - settle_dt) / G_DAYS_IN_YEAR

        coupon_times.append(flow_time)
        coupon_flows.append(cpn)

coupon_times = np.array(coupon_times)
coupon_flows = np.array(coupon_flows)


# ============================================================================
# 1.3 DISCOUNT FACTORS
# ============================================================================

y = 0.05

times = np.linspace(
    0.0,
    10.0,
    21,
)

dfs = np.power(
    1.0 + y / 2.0,
    -times * 2.0,
)


# ============================================================================
# 1.4 HULL-WHITE MODEL
# ============================================================================

sigma = 0.0125
a = 0.10

model = HWTree(
    sigma,
    a,
    None,
)

t_exp = (expiry_dt - settle_dt) / G_DAYS_IN_YEAR

t_mat = (maturity_dt - settle_dt) / G_DAYS_IN_YEAR


# ============================================================================
# 1.5 JAMSHIDIAN DECOMPOSITION
# ============================================================================
#
# The Jamshidian method returns:
#
#     (call_value, put_value)
#
# Unpack the two values explicitly.
# ============================================================================

value_jamshidian = model.european_bond_option_jamshidian(
    t_exp,
    strike_price,
    face,
    coupon_times,
    coupon_flows,
    times,
    dfs,
)

jam_call, jam_put = value_jamshidian


# ============================================================================
# 1.6 TREE VALUATION
# ============================================================================
#
# Build the Hull-White tree and value the same European bond option.
#
# The tree calculation also returns:
#
#     (call_value, put_value)
# ============================================================================

num_time_steps = 100

model.num_time_steps = num_time_steps

model.build_tree(
    t_mat,
    times,
    dfs,
)

exercise_type = ExerciseTypes.EUROPEAN

value_tree = model.bond_option(
    t_exp,
    strike_price,
    face,
    coupon_times,
    coupon_flows,
    exercise_type,
)

tree_call, tree_put = value_tree


# ============================================================================
# 1.7 RESULTS
# ============================================================================

print("\n" + "-" * 72)

print(f"{'METHOD':<26}" f"{'CALL':>16}" f"{'PUT':>16}")

print("-" * 72)

print(f"{'Jamshidian':<26}" f"{jam_call:16.6f}" f"{jam_put:16.6f}")

print(f"{'Hull-White Tree':<26}" f"{tree_call:16.6f}" f"{tree_put:16.6f}")

print("-" * 72)

print(f"{'Tree - Jamshidian':<26}" f"{tree_call - jam_call:16.6f}" f"{tree_put - jam_put:16.6f}")

print("-" * 72)


# ============================================================================
# 2. BOND OPTION - NEAR-ZERO-VOLATILITY CONVERGENCE
# ============================================================================
#
# At very low volatility, option values should approach their deterministic
# values.
#
# European option values can therefore be compared with discounted
# intrinsic values calculated from the forward clean bond price.
# ============================================================================

print("\n" + "=" * 100)
print("2. BOND OPTION - NEAR-ZERO-VOLATILITY CONVERGENCE")
print("=" * 100)


# ============================================================================
# 2.1 BOND AND CURVE
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

expiry_dt = Date(1, 12, 2021)


# ============================================================================
# 2.2 REFERENCE VALUES
# ============================================================================

df_expiry = discount_curve.df(
    expiry_dt,
)

fwd_clean_value = bond.clean_price_from_discount_curve(
    expiry_dt,
    discount_curve,
)

spot_clean_value = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Spot clean bond price    : {spot_clean_value:.6f}")
print(f"Forward clean bond price : {fwd_clean_value:.6f}")
print(f"Expiry discount factor   : {df_expiry:.6f}")


# ============================================================================
# 2.3 CONVERGENCE TEST
# ============================================================================

num_time_steps = range(
    100,
    400,
    100,
)

strike_prices = [
    90.0,
    120.0,
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

        model = HWTree(
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
# 3. BOND OPTION VALUES BY STRIKE
# ============================================================================
#
# Value European and American calls and puts for several strikes using the
# Hull-White model.
# ============================================================================

print("\n" + "=" * 100)
print("3. BOND OPTION VALUES BY STRIKE")
print("=" * 100)


# ============================================================================
# 3.1 BOND AND CURVE
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

times = np.linspace(
    0.0,
    10.0,
    21,
)

dfs = np.exp(
    -0.05 * times,
)

dates = settle_dt.add_years(
    times,
)

discount_curve = DiscountCurve(
    settle_dt,
    dates,
    dfs,
)

expiry_dt = settle_dt.add_tenor("18m")

strikes = [
    80.0,
    90.0,
    100.0,
    110.0,
    120.0,
]

price = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Clean bond price : {price:.6f}")


# ============================================================================
# 3.2 OPTION VALUES
# ============================================================================

num_time_steps = 50

print("\n" + "-" * 88)

print(f"{'OPTION TYPE':<22}" f"{'STRIKE':>10}" f"{'SIGMA':>12}" f"{'MEAN REV':>12}" f"{'VALUE':>16}")

print("-" * 88)

option_specs = [
    (
        "EUROPEAN CALL",
        OptionTypes.EUROPEAN_CALL,
        0.01,
    ),
    (
        "AMERICAN CALL",
        OptionTypes.AMERICAN_CALL,
        0.01,
    ),
    (
        "EUROPEAN PUT",
        OptionTypes.EUROPEAN_PUT,
        0.01,
    ),
    (
        "AMERICAN PUT",
        OptionTypes.AMERICAN_PUT,
        0.02,
    ),
]

a = 0.10

for option_label, option_type, sigma in option_specs:

    for strike_price in strikes:

        bond_option = BondOption(
            bond,
            expiry_dt,
            strike_price,
            option_type,
        )

        model = HWTree(
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
# 4. EUROPEAN OPTION CALCULATION CONVERGENCE
# ============================================================================
#
# Compare alternative Hull-White approaches for European bond options.
#
# The expiry is deliberately placed between coupon dates.
#
# For puts, compare the standard calculation with EXPIRY_ONLY.
# For calls, compare the standard calculation with EXPIRY_TREE.
# ============================================================================

print("\n" + "=" * 100)
print("4. EUROPEAN OPTION CALCULATION CONVERGENCE")
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

issue_dt = Date(1, 12, 2015)
maturity_dt = Date(1, 12, 2020)

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

# Put expiry in the middle of a coupon period.

expiry_dt = Date(1, 3, 2020)

strike_price = 100.0

sigma = 0.05
a = 0.10

time_steps = range(
    100,
    400,
    100,
)


# ============================================================================
# 4.2 CONVERGENCE
# ============================================================================

print("\n" + "-" * 94)

print(
    f"{'STEPS':>8}" f"{'PUT STD':>16}" f"{'PUT EXPIRY':>16}" f"{'CALL STD':>16}" f"{'CALL EXPIRY':>16}" f"{'TIME':>16}"
)

print("-" * 94)

for num_time_steps in time_steps:

    start = time.perf_counter()

    euro_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_PUT,
    )

    model_put_standard = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    value_put_standard = euro_put.value(
        settle_dt,
        discount_curve,
        model_put_standard,
    )

    model_put_expiry = HWTree(
        sigma,
        a,
        num_time_steps,
        HWEuropeanCalcTypes.EXPIRY_ONLY,
    )

    value_put_expiry = euro_put.value(
        settle_dt,
        discount_curve,
        model_put_expiry,
    )

    euro_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    model_call_standard = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    value_call_standard = euro_call.value(
        settle_dt,
        discount_curve,
        model_call_standard,
    )

    model_call_expiry = HWTree(
        sigma,
        a,
        num_time_steps,
        HWEuropeanCalcTypes.EXPIRY_TREE,
    )

    value_call_expiry = euro_call.value(
        settle_dt,
        discount_curve,
        model_call_expiry,
    )

    elapsed = time.perf_counter() - start

    print(
        f"{num_time_steps:8d}"
        f"{value_put_standard:16.6f}"
        f"{value_put_expiry:16.6f}"
        f"{value_call_standard:16.6f}"
        f"{value_call_expiry:16.6f}"
        f"{elapsed:16.6f}"
    )


# ============================================================================
# 5. AMERICAN VERSUS EUROPEAN OPTION CONVERGENCE
# ============================================================================
#
# Compare American and European puts and calls as the number of Hull-White
# tree time steps is increased.
#
# European options use the specialised expiry calculations retained from
# the original example.
# ============================================================================

print("\n" + "=" * 100)
print("5. AMERICAN VERSUS EUROPEAN OPTION CONVERGENCE")
print("=" * 100)


# ============================================================================
# 5.1 BOND AND CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

discount_curve = FlatDiscountCurve(
    settle_dt,
    0.05,
)

issue_dt = Date(1, 9, 2014)
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

sigma = 0.05
a = 0.10

time_steps = range(
    100,
    400,
    100,
)


# ============================================================================
# 5.2 CONVERGENCE
# ============================================================================

print("\n" + "-" * 94)

print(f"{'STEPS':>8}" f"{'AMER PUT':>16}" f"{'EUR PUT':>16}" f"{'AMER CALL':>16}" f"{'EUR CALL':>16}" f"{'TIME':>16}")

print("-" * 94)

for num_time_steps in time_steps:

    start = time.perf_counter()

    amer_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_PUT,
    )

    model_amer_put = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    value_amer_put = amer_put.value(
        settle_dt,
        discount_curve,
        model_amer_put,
    )

    euro_put = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_PUT,
    )

    model_euro_put = HWTree(
        sigma,
        a,
        num_time_steps,
        HWEuropeanCalcTypes.EXPIRY_ONLY,
    )

    value_euro_put = euro_put.value(
        settle_dt,
        discount_curve,
        model_euro_put,
    )

    amer_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.AMERICAN_CALL,
    )

    model_amer_call = HWTree(
        sigma,
        a,
        num_time_steps,
    )

    value_amer_call = amer_call.value(
        settle_dt,
        discount_curve,
        model_amer_call,
    )

    euro_call = BondOption(
        bond,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
    )

    model_euro_call = HWTree(
        sigma,
        a,
        num_time_steps,
        HWEuropeanCalcTypes.EXPIRY_TREE,
    )

    value_euro_call = euro_call.value(
        settle_dt,
        discount_curve,
        model_euro_call,
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
# 6. DETAILED TREE CONVERGENCE
# ============================================================================
#
# Examine convergence of European and American call and put values using
# increasingly fine Hull-White trees.
#
# The convergence series are plotted automatically.
# ============================================================================

print("\n" + "=" * 100)
print("6. DETAILED TREE CONVERGENCE")
print("=" * 100)


# ============================================================================
# 6.1 BOND AND CURVE
# ============================================================================

settle_dt = Date(1, 12, 2019)

discount_curve = FlatDiscountCurve(
    settle_dt,
    0.05,
)

issue_dt = Date(1, 12, 2015)
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

expiry_dt = settle_dt.add_tenor("18m")

spot_value = bond.clean_price_from_discount_curve(
    settle_dt,
    discount_curve,
)

print(f"Clean bond price : {spot_value:.6f}")


# ============================================================================
# 6.2 OPTION AND MODEL INPUTS
# ============================================================================

strike_price = 102.0

sigma = 0.01
a = 0.10

print(f"Strike           : {strike_price:.2f}")
print(f"HW volatility    : {sigma:.6f}")
print(f"HW mean reversion: {a:.6f}")


# ============================================================================
# 6.3 TREE CONVERGENCE
# ============================================================================

num_steps_vector = range(
    100,
    400,
    100,
)

values_euro_call = []
values_amer_call = []
values_euro_put = []
values_amer_put = []

print("\n" + "-" * 92)

print(f"{'STEPS':>8}" f"{'EUR CALL':>16}" f"{'AMER CALL':>16}" f"{'EUR PUT':>16}" f"{'AMER PUT':>16}" f"{'TIME':>16}")

print("-" * 92)

for num_steps in num_steps_vector:

    model = HWTree(
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
# 6.4 CONVERGENCE PLOTS
# ============================================================================
#
# Plot all four option-value convergence series.
#
# The plots are always generated when the example is run.
# ============================================================================

steps = list(num_steps_vector)

plt.figure()
plt.plot(
    steps,
    values_euro_call,
    marker="o",
    label="European Call",
)
plt.title("European Call Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_call,
    marker="o",
    label="American Call",
)
plt.title("American Call Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_euro_put,
    marker="o",
    label="European Put",
)
plt.title("European Put Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.figure()
plt.plot(
    steps,
    values_amer_put,
    marker="o",
    label="American Put",
)
plt.title("American Put Convergence - Hull-White Model")
plt.xlabel("Number of Hull-White Time Steps")
plt.ylabel("Option Value")
plt.legend()
plt.grid()

plt.show()
