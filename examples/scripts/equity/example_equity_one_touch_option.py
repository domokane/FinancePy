# ============================================================================
# FINANCEPY EXAMPLES - EquityOneTouchOption
# ============================================================================
#
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
#

import numpy as np
import matplotlib.pyplot as plt

from financepy.products.equity.equity_one_touch_option import EquityOneTouchOption
from financepy.utils.global_types import TouchOptionTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.utils.date import Date


value_dt = Date(1, 1, 2016)
expiry_dt = Date(2, 7, 2016)

interest_rate = 0.10
dividend_yield = 0.03
volatility = 0.20
barrier_level = 100.0
payment_size = 15.0

model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

num_paths = 10000
num_observations_per_year = 252

# ============================================================================
# 1. CASH ONE-TOUCH OPTIONS
# ============================================================================
# Down-barrier contracts start above the barrier; up-barrier contracts start
# below the barrier. Analytic values use continuous barrier observation, while
# the Monte Carlo calculation observes the barrier on its simulation grid.

print("\n" + "=" * 78)
print("1. CASH ONE-TOUCH OPTIONS")
print("=" * 78)

cash_cases = [
    (105.0, TouchOptionTypes.DOWN_AND_IN_CASH_AT_HIT),
    (105.0, TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY),
    (105.0, TouchOptionTypes.DOWN_AND_OUT_CASH_OR_NOTHING),
    (95.0, TouchOptionTypes.UP_AND_IN_CASH_AT_HIT),
    (95.0, TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY),
    (95.0, TouchOptionTypes.UP_AND_OUT_CASH_OR_NOTHING),
]

print(
    f"{'TYPE':>48s}"
    f"{'SPOT':>10s}"
    f"{'ANALYTIC':>14s}"
    f"{'MC':>14s}"
    f"{'MC-AN':>14s}"
)
print("-" * 100)

for stock_price, option_type in cash_cases:
    option = EquityOneTouchOption(
        expiry_dt, option_type, barrier_level, payment_size
    )

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_observations_per_year,
        num_paths,
    )

    print(
        f"{str(option_type):>48s}"
        f"{stock_price:10.4f}"
        f"{value:14.8f}"
        f"{value_mc:14.8f}"
        f"{value_mc - value:14.8f}"
    )


# ============================================================================
# 2. ASSET ONE-TOUCH OPTIONS
# ============================================================================

print("\n" + "=" * 78)
print("2. ASSET ONE-TOUCH OPTIONS")
print("=" * 78)

asset_cases = [
    (105.0, TouchOptionTypes.DOWN_AND_IN_ASSET_AT_HIT),
    (105.0, TouchOptionTypes.DOWN_AND_IN_ASSET_AT_EXPIRY),
    (105.0, TouchOptionTypes.DOWN_AND_OUT_ASSET_OR_NOTHING),
    (95.0, TouchOptionTypes.UP_AND_IN_ASSET_AT_HIT),
    (95.0, TouchOptionTypes.UP_AND_IN_ASSET_AT_EXPIRY),
    (95.0, TouchOptionTypes.UP_AND_OUT_ASSET_OR_NOTHING),
]

print(
    f"{'TYPE':>48s}"
    f"{'SPOT':>10s}"
    f"{'ANALYTIC':>14s}"
    f"{'MC':>14s}"
    f"{'MC-AN':>14s}"
)
print("-" * 100)

for stock_price, option_type in asset_cases:
    option = EquityOneTouchOption(
        expiry_dt, option_type, barrier_level
    )

    value = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_observations_per_year,
        num_paths,
    )

    print(
        f"{str(option_type):>48s}"
        f"{stock_price:10.4f}"
        f"{value:14.8f}"
        f"{value_mc:14.8f}"
        f"{value_mc - value:14.8f}"
    )


# ============================================================================
# 3. MONTE CARLO PATH CONVERGENCE
# ============================================================================
# Hold the observation grid fixed and increase only the number of paths.
# For each path count, run several independent replications using different
# seeds. This separates Monte Carlo sampling error from any persistent
# discrete-observation bias.

print("\n" + "=" * 78)
print("3. MONTE CARLO PATH CONVERGENCE")
print("=" * 78)

stock_price = 105.0
option_type = TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY

option = EquityOneTouchOption(
    expiry_dt, option_type, barrier_level, payment_size
)

analytic_value = option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

num_paths_list = [
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
]

num_reps = 20
base_seed = 4242

mc_means = []
mc_stds = []
mc_mean_ses = []
mc_biases = []

print(f"Option type       : {option_type}")
print(f"Analytic value    : {analytic_value:.10f}")
print(f"Observations/year : {num_observations_per_year}")
print(f"Replications      : {num_reps}")
print()
print(
    f"{'PATHS':>10s}"
    f"{'MC MEAN':>16s}"
    f"{'MC STD':>16s}"
    f"{'MEAN SE':>16s}"
    f"{'BIAS':>16s}"
    f"{'BIAS/SE':>12s}"
)
print("-" * 86)

for num_paths in num_paths_list:
    values = np.empty(num_reps)

    for rep in range(num_reps):
        values[rep] = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_observations_per_year,
            num_paths,
            seed=base_seed + rep,
        )

    mc_mean = np.mean(values)
    mc_std = np.std(values, ddof=1)
    mc_mean_se = mc_std / np.sqrt(num_reps)
    mc_bias = mc_mean - analytic_value
    bias_over_se = mc_bias / mc_mean_se if mc_mean_se > 0.0 else np.nan

    mc_means.append(mc_mean)
    mc_stds.append(mc_std)
    mc_mean_ses.append(mc_mean_se)
    mc_biases.append(mc_bias)

    print(
        f"{num_paths:10d}"
        f"{mc_mean:16.8f}"
        f"{mc_std:16.8f}"
        f"{mc_mean_se:16.8f}"
        f"{mc_bias:16.8f}"
        f"{bias_over_se:12.4f}"
    )

mc_means = np.asarray(mc_means)
mc_stds = np.asarray(mc_stds)
mc_mean_ses = np.asarray(mc_mean_ses)
mc_biases = np.asarray(mc_biases)

plt.figure(figsize=(10, 6))
plt.errorbar(
    num_paths_list,
    mc_means,
    yerr=mc_mean_ses,
    marker="o",
    capsize=4,
    label="Monte Carlo Mean +/- SE",
)
plt.axhline(
    analytic_value,
    linestyle="--",
    label="Analytic",
)
plt.xscale("log")
plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Option Value")
plt.title("One-Touch Option: Monte Carlo Path Convergence")
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.errorbar(
    num_paths_list,
    mc_biases,
    yerr=mc_mean_ses,
    marker="o",
    capsize=4,
    label="MC Mean - Analytic +/- SE",
)
plt.axhline(0.0, linestyle="--", label="Zero Error")
plt.xscale("log")
plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("MC Mean - Analytic Value")
plt.title("One-Touch Option: Monte Carlo Sampling Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. BARRIER OBSERVATION-FREQUENCY CONVERGENCE
# ============================================================================
# Hold the number of paths fixed and refine the simulation grid. A discrete
# grid can miss barrier crossings between observations, so this test is
# distinct from ordinary Monte Carlo path convergence.

print("\n" + "=" * 78)
print("4. BARRIER OBSERVATION-FREQUENCY CONVERGENCE")
print("=" * 78)

observations_per_year_list = [
    12,
    26,
    52,
    100,
    252,
    500,
    1000,
    2000,
    4000,
]

observation_num_paths = 100000
observation_values = []

print(f"Option type    : {option_type}")
print(f"Analytic value : {analytic_value:.10f}")
print(f"Paths          : {observation_num_paths}")
print()
print(
    f"{'OBS/YEAR':>12s}"
    f"{'MC VALUE':>16s}"
    f"{'ANALYTIC':>16s}"
    f"{'ERROR':>16s}"
)
print("-" * 60)

for observations_per_year in observations_per_year_list:
    value_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        observations_per_year,
        observation_num_paths,
    )

    observation_values.append(value_mc)

    print(
        f"{observations_per_year:12d}"
        f"{value_mc:16.8f}"
        f"{analytic_value:16.8f}"
        f"{value_mc - analytic_value:16.8f}"
    )

plt.figure(figsize=(10, 6))
plt.plot(
    observations_per_year_list,
    observation_values,
    marker="o",
    label="Monte Carlo",
)
plt.axhline(
    analytic_value,
    linestyle="--",
    label="Analytic",
)
plt.xlabel("Observations per Year")
plt.ylabel("Option Value")
plt.title("One-Touch Option: Observation-Frequency Convergence")
plt.grid(True)
plt.legend()
plt.show()

observation_errors = np.asarray(observation_values) - analytic_value

plt.figure(figsize=(10, 6))
plt.plot(
    observations_per_year_list,
    observation_errors,
    marker="o",
    label="MC - Analytic",
)
plt.axhline(0.0, linestyle="--", label="Zero Error")
plt.xlabel("Observations per Year")
plt.ylabel("MC Value - Analytic Value")
plt.title("One-Touch Option: Discrete-Observation Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. DOWN AND UP BARRIER OBSERVATION COMPARISON
# ============================================================================
# Compare representative down- and up-barrier contracts using the same
# discrete observation frequencies.

print("\n" + "=" * 78)
print("5. DOWN AND UP BARRIER OBSERVATIONS COMPARISON")
print("=" * 78)

comparison_cases = [
    (
        "DOWN-AND-IN CASH AT EXPIRY",
        105.0,
        TouchOptionTypes.DOWN_AND_IN_CASH_AT_EXPIRY,
    ),
    (
        "UP-AND-IN CASH AT EXPIRY",
        95.0,
        TouchOptionTypes.UP_AND_IN_CASH_AT_EXPIRY,
    ),
]

plt.figure(figsize=(10, 6))

for label, spot, touch_type in comparison_cases:
    comparison_option = EquityOneTouchOption(
        expiry_dt,
        touch_type,
        barrier_level,
        payment_size,
    )

    comparison_analytic = comparison_option.value(
        value_dt,
        spot,
        discount_curve,
        dividend_curve,
        model,
    )

    comparison_values = []

    for observations_per_year in observations_per_year_list:
        comparison_mc = comparison_option.value_mc(
            value_dt,
            spot,
            discount_curve,
            dividend_curve,
            model,
            observations_per_year,
            observation_num_paths,
        )

        comparison_values.append(comparison_mc)

    comparison_errors = (
        np.asarray(comparison_values) - comparison_analytic
    )

    print()
    print(label)
    print(f"Analytic value : {comparison_analytic:.10f}")
    print(
        f"{'OBS/YEAR':>12s}"
        f"{'MC VALUE':>16s}"
        f"{'ERROR':>16s}"
    )
    print("-" * 44)

    for observations_per_year, comparison_mc, error in zip(
        observations_per_year_list,
        comparison_values,
        comparison_errors,
    ):
        print(
            f"{observations_per_year:12d}"
            f"{comparison_mc:16.8f}"
            f"{error:16.8f}"
        )

    plt.plot(
        observations_per_year_list,
        comparison_errors,
        marker="o",
        label=label,
    )

plt.axhline(0.0, linestyle="--", label="Zero Error")
plt.xlabel("Observations per Year")
plt.ylabel("MC Value - Analytic Value")
plt.title("One-Touch Options: Barrier Observations Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("6. SUMMARY")
print("=" * 78)

print(f"Convergence option       : {option_type}")
print(f"Analytic value           : {analytic_value:.10f}")
print(f"Highest-path MC mean     : {mc_means[-1]:.10f}")
print(f"Highest-path bias        : {mc_biases[-1]:.10f}")
print(f"Finest-observation MC value: {observation_values[-1]:.10f}")
print(f"Finest-observation error   : {observation_errors[-1]:.10f}")


for seed in range(4242, 4252):

    v = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_observations_per_year,
        100000,
        seed=seed,
    )

    print(seed, v)
