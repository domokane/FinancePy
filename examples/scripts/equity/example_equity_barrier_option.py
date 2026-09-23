# ============================================================================
# FINANCEPY EXAMPLES - EquityBarrierOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

from financepy.models.black_scholes import BlackScholes

from financepy.products.equity.equity_barrier_option import BarrierTypes
from financepy.products.equity.equity_barrier_option import EquityBarrierOption
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


LINE = "=" * 78


# ============================================================================
# 1. EQUITY BARRIER OPTION
# ============================================================================
# What this section demonstrates:
#
# Compares the analytic barrier-option valuation with Monte Carlo valuation
# across all barrier types.
#
# Barrier options depend not only on the terminal stock price but also on
# whether the stock price has crossed a specified barrier during the life of
# the option.
#
# The comparison provides a useful numerical check that the Monte Carlo
# implementation is consistent with the analytic valuation.
# ============================================================================

print("\n" + LINE)
print("1. EQUITY BARRIER OPTION")
print(LINE)

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

model = BlackScholes(
    volatility,
)

num_obs_per_year = 252
num_paths = 10000
seed = 42

print(
    f"{'TYPE':>25}"
    f"{'K':>8}"
    f"{'B':>8}"
    f"{'S':>8}"
    f"{'VALUE':>14}"
    f"{'MC VALUE':>14}"
    f"{'DIFF':>14}"
    f"{'TIME':>10}"
)

print("-" * 101)

for opt_type in BarrierTypes:

    for stock_price in [80.0, 100.0, 120.0]:

        barrier = 110.0
        strike = 100.0

        option = EquityBarrierOption(
            expiry_dt,
            strike,
            opt_type,
            barrier,
            num_obs_per_year,
        )

        analytic_value = option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        start = time.time()

        mc_value = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_obs_per_year=num_obs_per_year,
            num_paths=num_paths,
            seed=seed,
        )

        end = time.time()

        elapsed = end - start
        difference = mc_value - analytic_value

        print(
            f"{str(opt_type):>25}"
            f"{strike:8.2f}"
            f"{barrier:8.2f}"
            f"{stock_price:8.2f}"
            f"{analytic_value:14.6f}"
            f"{mc_value:14.6f}"
            f"{difference:14.6f}"
            f"{elapsed:10.4f}"
        )

    for stock_price in [80.0, 100.0, 120.0]:

        strike = 110.0
        barrier = 100.0

        option = EquityBarrierOption(
            expiry_dt,
            strike,
            opt_type,
            barrier,
            num_obs_per_year,
        )

        analytic_value = option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        start = time.time()

        mc_value = option.value_mc(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_obs_per_year=num_obs_per_year,
            num_paths=num_paths,
            seed=seed,
        )

        end = time.time()

        elapsed = end - start
        difference = mc_value - analytic_value

        print(
            f"{str(opt_type):>25}"
            f"{strike:8.2f}"
            f"{barrier:8.2f}"
            f"{stock_price:8.2f}"
            f"{analytic_value:14.6f}"
            f"{mc_value:14.6f}"
            f"{difference:14.6f}"
            f"{elapsed:10.4f}"
        )


# ============================================================================
# 2. BARRIER OPTION GREEKS
# ============================================================================
# What this section demonstrates:
#
# Calculates the option value and selected Greeks for every barrier type.
#
# Barrier-option Greeks can change sharply when the stock price approaches
# the barrier because a relatively small movement in the underlying can
# materially change the probability of the barrier being touched.
# ============================================================================

print("\n" + LINE)
print("2. BARRIER OPTION GREEKS")
print(LINE)

stock_prices = [
    80.0,
    100.0,
    120.0,
]

strike = 100.0
barrier = 105.0

print(
    f"{'TYPE':>25}"
    f"{'S':>10}"
    f"{'VALUE':>14}"
    f"{'DELTA':>14}"
    f"{'VEGA':>14}"
    f"{'THETA':>14}"
)

print("-" * 91)

for opt_type in BarrierTypes:

    for stock_price in stock_prices:

        barrier_option = EquityBarrierOption(
            expiry_dt,
            strike,
            opt_type,
            barrier,
            num_obs_per_year,
        )

        option_value = barrier_option.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        delta = barrier_option.delta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        vega = barrier_option.vega(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        theta = barrier_option.theta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        print(
            f"{str(opt_type):>25}"
            f"{stock_price:10.2f}"
            f"{option_value:14.6f}"
            f"{delta:14.6f}"
            f"{vega:14.6f}"
            f"{theta:14.6f}"
        )


# ============================================================================
# 3. DOWN-AND-IN AND DOWN-AND-OUT CALL VALUES
# ============================================================================
# What this section demonstrates:
#
# A down-and-out call disappears if the stock touches the lower barrier.
#
# A down-and-in call becomes active if the same barrier is touched.
#
# The two contracts therefore divide the possible paths of the corresponding
# vanilla call into mutually exclusive sets. Their values should consequently
# add to approximately the vanilla call value.
#
# The plot also makes the behaviour close to the barrier visible.
# ============================================================================

print("\n" + LINE)
print("3. DOWN-AND-IN AND DOWN-AND-OUT CALL VALUES")
print(LINE)

strike = 100.0
barrier = 90.0

stock_grid = np.linspace(
    70.0,
    130.0,
    61,
)

down_in_values = []
down_out_values = []

down_in_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_IN_CALL,
    barrier,
    num_obs_per_year,
)

down_out_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

for stock_price in stock_grid:

    down_in_value = down_in_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    down_out_value = down_out_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    down_in_values.append(
        down_in_value,
    )

    down_out_values.append(
        down_out_value,
    )

down_in_values = np.asarray(
    down_in_values,
)

down_out_values = np.asarray(
    down_out_values,
)

plt.figure()

plt.plot(
    stock_grid,
    down_in_values,
    label="Down-and-In Call",
)

plt.plot(
    stock_grid,
    down_out_values,
    label="Down-and-Out Call",
)

plt.axvline(
    barrier,
    linestyle="--",
    linewidth=1.0,
    label="Barrier",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Down Barrier Call Value versus Stock Price"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. UP-AND-IN AND UP-AND-OUT CALL VALUES
# ============================================================================
# What this section demonstrates:
#
# The same complementary behaviour can be seen for an upper barrier.
#
# An up-and-out call loses value as the probability of hitting the upper
# barrier increases. The corresponding up-and-in call gains value.
# ============================================================================

print("\n" + LINE)
print("4. UP-AND-IN AND UP-AND-OUT CALL VALUES")
print(LINE)

strike = 100.0
barrier = 120.0

stock_grid = np.linspace(
    70.0,
    130.0,
    61,
)

up_in_values = []
up_out_values = []

up_in_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.UP_AND_IN_CALL,
    barrier,
    num_obs_per_year,
)

up_out_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.UP_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

for stock_price in stock_grid:

    up_in_value = up_in_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    up_out_value = up_out_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    up_in_values.append(
        up_in_value,
    )

    up_out_values.append(
        up_out_value,
    )

up_in_values = np.asarray(
    up_in_values,
)

up_out_values = np.asarray(
    up_out_values,
)

plt.figure()

plt.plot(
    stock_grid,
    up_in_values,
    label="Up-and-In Call",
)

plt.plot(
    stock_grid,
    up_out_values,
    label="Up-and-Out Call",
)

plt.axvline(
    barrier,
    linestyle="--",
    linewidth=1.0,
    label="Barrier",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Up Barrier Call Value versus Stock Price"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. KNOCK-IN / KNOCK-OUT PARITY
# ============================================================================
# What this section demonstrates:
#
# For otherwise identical contracts:
#
#       knock-in option + knock-out option = vanilla option
#
# because every stock-price path must either touch the barrier or not touch
# the barrier.
#
# This provides a particularly useful internal consistency check for the
# barrier-option implementation.
# ============================================================================

print("\n" + LINE)
print("5. KNOCK-IN / KNOCK-OUT PARITY")
print(LINE)

strike = 100.0
barrier = 90.0

vanilla_call = EquityVanillaOption(
    expiry_dt,
    strike,
    OptionTypes.EUROPEAN_CALL,
)

down_in_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_IN_CALL,
    barrier,
    num_obs_per_year,
)

down_out_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

print(
    f"{'S':>10}"
    f"{'VANILLA':>14}"
    f"{'DOWN-IN':>14}"
    f"{'DOWN-OUT':>14}"
    f"{'IN + OUT':>14}"
    f"{'ERROR':>14}"
)

print("-" * 80)

for stock_price in [
    95.0,
    100.0,
    105.0,
    110.0,
    120.0,
]:

    vanilla_value = vanilla_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    knock_in_value = down_in_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    knock_out_value = down_out_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    combined_value = (
        knock_in_value
        + knock_out_value
    )

    parity_error = (
        combined_value
        - vanilla_value
    )

    print(
        f"{stock_price:10.2f}"
        f"{vanilla_value:14.6f}"
        f"{knock_in_value:14.6f}"
        f"{knock_out_value:14.6f}"
        f"{combined_value:14.6f}"
        f"{parity_error:14.8f}"
    )


# ============================================================================
# 6. KNOCK-IN / KNOCK-OUT PARITY VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
#
# Extends the parity check in Section 5 across a continuous range of stock
# prices.
#
# For otherwise identical contracts:
#
#       knock-in option + knock-out option = vanilla option
#
# The parity error should therefore remain close to zero across the stock-price
# grid. Plotting the error provides a useful visual regression test for the
# barrier-option implementation.
# ============================================================================

print("\n" + LINE)
print("6. KNOCK-IN / KNOCK-OUT PARITY VERSUS STOCK PRICE")
print(LINE)

strike = 100.0
barrier = 90.0

stock_grid = np.linspace(
    90.5,
    130.0,
    80,
)

vanilla_call = EquityVanillaOption(
    expiry_dt,
    strike,
    OptionTypes.EUROPEAN_CALL,
)

down_in_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_IN_CALL,
    barrier,
    num_obs_per_year,
)

down_out_call = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

parity_errors = []

for stock_price in stock_grid:

    vanilla_value = vanilla_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    knock_in_value = down_in_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    knock_out_value = down_out_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    parity_error = (
        knock_in_value
        + knock_out_value
        - vanilla_value
    )

    parity_errors.append(
        parity_error,
    )

parity_errors = np.asarray(
    parity_errors,
)

print(
    f"{'Maximum Absolute Parity Error':<35}"
    f"{np.max(np.abs(parity_errors)):20.10f}"
)

plt.figure()

plt.plot(
    stock_grid,
    parity_errors,
)

plt.axhline(
    0.0,
    linestyle="--",
    linewidth=1.0,
)

plt.axvline(
    barrier,
    linestyle="--",
    linewidth=1.0,
    label="Barrier",
)

plt.xlabel("Stock Price")
plt.ylabel("Knock-In + Knock-Out - Vanilla")

plt.title(
    "Barrier Option Knock-In / Knock-Out Parity Error"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. DOWN-AND-OUT CALL VALUE VERSUS BARRIER
# ============================================================================
# What this section demonstrates:
#
# The barrier level is one of the defining risk factors of a barrier option.
#
# For a down-and-out call, raising the barrier towards the current stock price
# makes it increasingly likely that the option will be knocked out.
#
# Its value should therefore generally fall as the lower barrier rises toward
# the current stock price.
# ============================================================================

print("\n" + LINE)
print("7. DOWN-AND-OUT CALL VALUE VERSUS BARRIER")
print(LINE)

stock_price = 100.0
strike = 100.0

barrier_grid = np.linspace(
    60.0,
    99.0,
    40,
)

barrier_values = []

for barrier in barrier_grid:

    barrier_option = EquityBarrierOption(
        expiry_dt,
        strike,
        BarrierTypes.DOWN_AND_OUT_CALL,
        barrier,
        num_obs_per_year,
    )

    option_value = barrier_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    barrier_values.append(
        option_value,
    )

barrier_values = np.asarray(
    barrier_values,
)

plt.figure()

plt.plot(
    barrier_grid,
    barrier_values,
    marker="o",
)

plt.xlabel("Barrier Level")
plt.ylabel("Option Value")

plt.title(
    "Down-and-Out Call Value versus Barrier"
)

plt.grid(True)
plt.show()


# ============================================================================
# 8. DELTA NEAR THE BARRIER
# ============================================================================
# What this section demonstrates:
#
# Barrier options can have very different hedge behaviour from vanilla
# options near the barrier.
#
# This plot shows the delta of a down-and-out call as the underlying stock
# price approaches and moves away from the barrier.
# ============================================================================

print("\n" + LINE)
print("8. DELTA NEAR THE BARRIER")
print(LINE)

strike = 100.0
barrier = 90.0

stock_grid = np.linspace(
    90.5,
    120.0,
    60,
)

delta_values = []

barrier_option = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

for stock_price in stock_grid:

    delta = barrier_option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    delta_values.append(
        delta,
    )

delta_values = np.asarray(
    delta_values,
)

plt.figure()

plt.plot(
    stock_grid,
    delta_values,
)

plt.axvline(
    barrier,
    linestyle="--",
    linewidth=1.0,
    label="Barrier",
)

plt.xlabel("Stock Price")
plt.ylabel("Delta")

plt.title(
    "Down-and-Out Call Delta Near the Barrier"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. ANALYTIC VERSUS MONTE CARLO
# ============================================================================
# What this section demonstrates:
#
# Compares the analytic and Monte Carlo values of a single barrier contract
# over a range of stock prices.
#
# The two curves should be close. Differences arise from Monte Carlo sampling
# error and from the treatment of barrier monitoring.
# ============================================================================

print("\n" + LINE)
print("9. ANALYTIC VERSUS MONTE CARLO")
print(LINE)

strike = 100.0
barrier = 90.0

stock_grid = np.linspace(
    92.0,
    125.0,
    12,
)

analytic_values = []
mc_values = []

barrier_option = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

num_paths = 10000

print(
    f"{'S':>10}"
    f"{'ANALYTIC':>16}"
    f"{'MONTE CARLO':>16}"
    f"{'DIFFERENCE':>16}"
)

print("-" * 58)

for stock_price in stock_grid:

    analytic_value = barrier_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    mc_value = barrier_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_obs_per_year=num_obs_per_year,
        num_paths=num_paths,
        seed=seed,
    )

    analytic_values.append(
        analytic_value,
    )

    mc_values.append(
        mc_value,
    )

    print(
        f"{stock_price:10.2f}"
        f"{analytic_value:16.6f}"
        f"{mc_value:16.6f}"
        f"{mc_value - analytic_value:16.6f}"
    )

analytic_values = np.asarray(
    analytic_values,
)

mc_values = np.asarray(
    mc_values,
)

plt.figure()

plt.plot(
    stock_grid,
    analytic_values,
    marker="o",
    label="Analytic",
)

plt.plot(
    stock_grid,
    mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Barrier Option Analytic and Monte Carlo Values"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 10. MONTE CARLO ERROR VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
#
# When the analytic and Monte Carlo curves lie almost on top of one another,
# plotting their difference makes the numerical error much easier to see.
# ============================================================================

print("\n" + LINE)
print("10. MONTE CARLO ERROR VERSUS STOCK PRICE")
print(LINE)

mc_errors = (
    mc_values
    - analytic_values
)

plt.figure()

plt.plot(
    stock_grid,
    mc_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Monte Carlo - Analytic")

plt.title(
    "Barrier Option Monte Carlo Valuation Error"
)

plt.grid(True)
plt.show()


# ============================================================================
# 11. MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
#
# Monte Carlo valuation contains sampling error. Increasing the number of
# simulated paths should generally cause the estimate to become more stable
# around the analytic value.
#
# Using the same random seed makes this experiment reproducible.
# ============================================================================

print("\n" + LINE)
print("11. MONTE CARLO CONVERGENCE")
print(LINE)

stock_price = 100.0
strike = 100.0
barrier = 90.0

barrier_option = EquityBarrierOption(
    expiry_dt,
    strike,
    BarrierTypes.DOWN_AND_OUT_CALL,
    barrier,
    num_obs_per_year,
)

analytic_value = barrier_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

path_counts = np.array(
    [
        100,
        250,
        500,
        1000,
        2500,
        5000,
        10000,
        25000,
    ]
)

mc_convergence_values = []

print(
    f"{'PATHS':>12}"
    f"{'MC VALUE':>16}"
    f"{'ANALYTIC':>16}"
    f"{'ERROR':>16}"
)

print("-" * 60)

for num_paths in path_counts:

    mc_value = barrier_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_obs_per_year=num_obs_per_year,
        num_paths=int(num_paths),
        seed=seed,
    )

    mc_convergence_values.append(
        mc_value,
    )

    print(
        f"{num_paths:12d}"
        f"{mc_value:16.6f}"
        f"{analytic_value:16.6f}"
        f"{mc_value - analytic_value:16.6f}"
    )

mc_convergence_values = np.asarray(
    mc_convergence_values,
)

plt.figure()

plt.plot(
    path_counts,
    mc_convergence_values,
    marker="o",
    label="Monte Carlo",
)

plt.axhline(
    analytic_value,
    linestyle="--",
    linewidth=1.0,
    label="Analytic",
)

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Option Value")

plt.title(
    "Barrier Option Monte Carlo Convergence"
)

plt.grid(True)
plt.legend()
plt.show()
