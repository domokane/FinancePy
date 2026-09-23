# ============================================================================
# FINANCEPY EXAMPLES - EquityAsianOption
# ============================================================================

import time
import numpy as np
import matplotlib.pyplot as plt


from financepy.utils.date import Date

from financepy.models.black_scholes import BlackScholes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_asian_option import (
    AsianOptionValuationTypes,
    EquityAsianOption,
)
from financepy.utils.global_types import OptionTypes


# ============================================================================
# 1. ASIAN OPTION VALUATION METHODS
# ============================================================================
# What this section demonstrates:
# Values an arithmetic-average Asian call using the different valuation
# methods available in FinancePy.
#
# The geometric-average option has an analytic solution and is useful as a
# reference calculation. Turnbull-Wakeman and Curran are approximations for
# arithmetic-average options, while the Monte Carlo methods simulate the
# underlying equity paths directly.
#
# Comparing the methods illustrates both model approximation error and
# Monte Carlo sampling error.
# ============================================================================

print("\n" + "=" * 78)
print("1. ASIAN OPTION VALUATION METHODS")
print("=" * 78)

value_dt = Date(1, 1, 2014)
start_averaging_dt = Date(1, 6, 2014)
expiry_dt = Date(1, 1, 2015)

stock_price = 100.0
strike_price = 100.0

volatility = 0.20
interest_rate = 0.30
dividend_yield = 0.10

num_obs_per_year = 120
num_paths = 1000
seed = 1991

accrued_average = stock_price * 1.10

model = BlackScholes(volatility)

print(f"{'METHOD':<25s}{'VALUE':>15s}")
print("-" * 40)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

asian_option = EquityAsianOption(
    start_averaging_dt,
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    num_obs_per_year,
)

value_geometric = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.GEOMETRIC,
    accrued_average,
)

print(f"{'Geometric':<25s}{value_geometric:15.8f}")

value_turnbull = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    accrued_average,
)

print(f"{'Turnbull-Wakeman':<25s}{value_turnbull:15.8f}")

value_curran = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.CURRAN,
    accrued_average,
)

print(f"{'Curran':<25s}{value_curran:15.8f}")

value_mc_fast = asian_option.value_mc_fast(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

print(f"{'Monte Carlo Fast':<25s}{value_mc_fast:15.8f}")

value_mc = asian_option.value_mc(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

print(f"{'Monte Carlo':<25s}{value_mc:15.8f}")


# ============================================================================
# 2. COMPARE ARITHMETIC ASIAN APPROXIMATIONS
# ============================================================================
# What this section demonstrates:
# Compares the Turnbull-Wakeman and Curran approximations with a Monte Carlo
# valuation.
#
# Arithmetic-average Asian options do not generally have the same simple
# closed-form solution as geometric-average Asian options. Approximation
# methods therefore provide a fast alternative to simulation.
#
# The difference relative to Monte Carlo gives an indication of the size of
# the approximation error for this particular set of market inputs.
# ============================================================================

print("\n" + "=" * 78)
print("2. COMPARE ARITHMETIC ASIAN APPROXIMATIONS")
print("=" * 78)

print(
    f"{'METHOD':<25s}"
    f"{'VALUE':>15s}"
    f"{'DIFF VS MC':>18s}"
)

print("-" * 58)

print(
    f"{'Turnbull-Wakeman':<25s}"
    f"{value_turnbull:15.8f}"
    f"{value_turnbull - value_mc:18.8f}"
)

print(
    f"{'Curran':<25s}"
    f"{value_curran:15.8f}"
    f"{value_curran - value_mc:18.8f}"
)

print(
    f"{'Monte Carlo':<25s}"
    f"{value_mc:15.8f}"
    f"{0.0:18.8f}"
)


# ============================================================================
# 3. MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Examines how the Monte Carlo Asian option value changes as the number of
# simulated paths increases.
#
# Monte Carlo valuation contains sampling error. Increasing the number of
# paths should reduce that error, although convergence is not necessarily
# monotonic.
#
# The analytic approximation values are also shown as fixed reference levels.
# ============================================================================

print("\n" + "=" * 78)
print("3. MONTE CARLO CONVERGENCE")
print("=" * 78)

num_paths_list = [
    1000,
    2000,
    5000,
    6000,
    8000,
    10000,
]

mc_values = []
mc_fast_values = []

print(
    f"{'PATHS':>10s}"
    f"{'MC':>18s}"
    f"{'MC FAST':>18s}"
    f"{'TW':>18s}"
    f"{'CURRAN':>18s}"
)

print("-" * 82)

for paths in num_paths_list:

    value_mc_test = asian_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    value_mc_fast_test = asian_option.value_mc_fast(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    mc_values.append(value_mc_test)
    mc_fast_values.append(value_mc_fast_test)

    print(
        f"{paths:10d}"
        f"{value_mc_test:18.8f}"
        f"{value_mc_fast_test:18.8f}"
        f"{value_turnbull:18.8f}"
        f"{value_curran:18.8f}"
    )


# ============================================================================
# 4. PLOT MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Plots the Monte Carlo estimates against the number of simulation paths.
#
# The horizontal Turnbull-Wakeman and Curran lines provide deterministic
# reference values. The simulation estimates should become progressively
# more stable as the number of paths increases.
# ============================================================================

print("\n" + "=" * 78)
print("4. PLOT MONTE CARLO CONVERGENCE")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    num_paths_list,
    mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.plot(
    num_paths_list,
    mc_fast_values,
    marker="o",
    label="Monte Carlo Fast",
)

plt.axhline(
    value_turnbull,
    linestyle="--",
    label="Turnbull-Wakeman",
)

plt.axhline(
    value_curran,
    linestyle="--",
    label="Curran",
)

plt.xlabel("Number of Paths")
plt.ylabel("Asian Option Value")
plt.title("Asian Option Monte Carlo Convergence")

plt.xscale("log")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. MONTE CARLO TIMINGS
# ============================================================================
# What this section demonstrates:
# Measures the computational cost of the Monte Carlo implementations as the
# number of simulated paths increases.
#
# Monte Carlo accuracy comes at a computational cost. This comparison shows
# how execution time grows as additional paths are used.
# ============================================================================

print("\n" + "=" * 78)
print("5. MONTE CARLO TIMINGS")
print("=" * 78)

mc_times = []
mc_fast_times = []

print(
    f"{'PATHS':>10s}"
    f"{'MC VALUE':>18s}"
    f"{'MC TIME':>15s}"
    f"{'FAST VALUE':>18s}"
    f"{'FAST TIME':>15s}"
)

print("-" * 76)

for paths in num_paths_list:

    start = time.time()

    value_mc_test = asian_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    end = time.time()

    mc_time = end - start

    start = time.time()

    value_mc_fast_test = asian_option.value_mc_fast(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    end = time.time()

    mc_fast_time = end - start

    mc_times.append(mc_time)
    mc_fast_times.append(mc_fast_time)

    print(
        f"{paths:10d}"
        f"{value_mc_test:18.8f}"
        f"{mc_time:15.6f}"
        f"{value_mc_fast_test:18.8f}"
        f"{mc_fast_time:15.6f}"
    )


# ============================================================================
# 6. PLOT MONTE CARLO TIMINGS
# ============================================================================
# What this section demonstrates:
# Shows the relationship between the number of simulated paths and execution
# time.
#
# This complements the convergence plot by illustrating the trade-off between
# numerical stability and computational cost.
# ============================================================================

print("\n" + "=" * 78)
print("6. PLOT MONTE CARLO TIMINGS")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    num_paths_list,
    mc_times,
    marker="o",
    label="Monte Carlo",
)

plt.plot(
    num_paths_list,
    mc_fast_times,
    marker="o",
    label="Monte Carlo Fast",
)

plt.xlabel("Number of Paths")
plt.ylabel("Calculation Time (seconds)")
plt.title("Asian Option Monte Carlo Calculation Time")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. OPTION VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how the Asian call value changes as the current stock price changes.
#
# A call option should generally become more valuable as the underlying stock
# price increases. Comparing the approximation and Monte Carlo methods also
# shows whether their agreement changes across moneyness.
# ============================================================================

print("\n" + "=" * 78)
print("7. OPTION VALUE VERSUS STOCK PRICE")
print("=" * 78)

stock_prices = np.linspace(
    70.0,
    130.0,
    13,
)

turnbull_stock_values = []
curran_stock_values = []
mc_stock_values = []

num_paths_stock = 10000

print(
    f"{'STOCK':>10s}"
    f"{'TURNBULL':>18s}"
    f"{'CURRAN':>18s}"
    f"{'MONTE CARLO':>18s}"
)

print("-" * 64)

for stock in stock_prices:

    accrued_average_stock = stock * 1.10

    turnbull_value = asian_option.value(
        value_dt,
        stock,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_average_stock,
    )

    curran_value = asian_option.value(
        value_dt,
        stock,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.CURRAN,
        accrued_average_stock,
    )

    mc_value = asian_option.value_mc(
        value_dt,
        stock,
        discount_curve,
        dividend_curve,
        model,
        num_paths_stock,
        seed,
        accrued_average_stock,
    )

    turnbull_stock_values.append(turnbull_value)
    curran_stock_values.append(curran_value)
    mc_stock_values.append(mc_value)

    print(
        f"{stock:10.2f}"
        f"{turnbull_value:18.8f}"
        f"{curran_value:18.8f}"
        f"{mc_value:18.8f}"
    )

turnbull_stock_values = np.array(turnbull_stock_values)
curran_stock_values = np.array(curran_stock_values)
mc_stock_values = np.array(mc_stock_values)

# ============================================================================
# 8. PLOT OPTION VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Visualises the Asian call's exposure to the underlying equity price and
# compares the different valuation methods across moneyness.
# ============================================================================

print("\n" + "=" * 78)
print("8. PLOT OPTION VALUE VERSUS STOCK PRICE")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    turnbull_stock_values,
    marker="o",
    label="Turnbull-Wakeman",
)

plt.plot(
    stock_prices,
    curran_stock_values,
    marker="o",
    label="Curran",
)

plt.plot(
    stock_prices,
    mc_stock_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Stock Price")
plt.ylabel("Asian Call Value")
plt.title("Asian Option Value versus Stock Price")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# PLOT DIFFERENCES BETWEEN VALUATION METHODS
# ============================================================================
#
# The three valuation methods are very close, so their curves overlap in the
# main plot. Plotting the differences makes the approximation errors visible.
# ============================================================================

plt.figure()

plt.plot(
    stock_prices,
    curran_stock_values - turnbull_stock_values,
    marker="o",
    label="Curran - Turnbull-Wakeman",
)

plt.plot(
    stock_prices,
    mc_stock_values - turnbull_stock_values,
    marker="o",
    label="Monte Carlo - Turnbull-Wakeman",
)

plt.axhline(
    0.0,
    linestyle="--",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Value Difference")
plt.title("Asian Option Valuation Method Differences")

plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 9. OPTION VALUE VERSUS VOLATILITY
# ============================================================================
# What this section demonstrates:
# Shows how the Asian option value responds to changes in equity volatility.
#
# Greater volatility generally increases option value because the holder
# benefits from favourable upside outcomes while downside exposure is limited
# by the option payoff.
#
# Asian options are less sensitive to volatility than otherwise comparable
# vanilla options because averaging reduces the variability of the effective
# underlying price.
# ============================================================================

print("\n" + "=" * 78)
print("9. OPTION VALUE VERSUS VOLATILITY")
print("=" * 78)

volatilities = np.linspace(
    0.05,
    0.60,
    12,
)

turnbull_vol_values = []
curran_vol_values = []
mc_vol_values = []

print(
    f"{'VOLATILITY':>12s}"
    f"{'TURNBULL':>18s}"
    f"{'CURRAN':>18s}"
    f"{'MONTE CARLO':>18s}"
)

print("-" * 66)

for vol in volatilities:

    vol_model = BlackScholes(vol)

    turnbull_value = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        vol_model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_average,
    )

    curran_value = asian_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        vol_model,
        AsianOptionValuationTypes.CURRAN,
        accrued_average,
    )

    mc_value = asian_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        vol_model,
        10000,
        seed,
        accrued_average,
    )

    turnbull_vol_values.append(turnbull_value)
    curran_vol_values.append(curran_value)
    mc_vol_values.append(mc_value)

    print(
        f"{vol:12.4f}"
        f"{turnbull_value:18.8f}"
        f"{curran_value:18.8f}"
        f"{mc_value:18.8f}"
    )


# ============================================================================
# 10. PLOT OPTION VALUE VERSUS VOLATILITY
# ============================================================================
# What this section demonstrates:
# Visualises the volatility sensitivity of the Asian option and shows how
# closely the analytic approximations track the Monte Carlo valuation.
# ============================================================================

print("\n" + "=" * 78)
print("10. PLOT OPTION VALUE VERSUS VOLATILITY")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    volatilities,
    turnbull_vol_values,
    marker="o",
    label="Turnbull-Wakeman",
)

plt.plot(
    volatilities,
    curran_vol_values,
    marker="o",
    label="Curran",
)

plt.plot(
    volatilities,
    mc_vol_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Volatility")
plt.ylabel("Asian Call Value")
plt.title("Asian Option Value versus Volatility")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 11. TIME EVOLUTION
# ============================================================================
# What this section demonstrates:
# Values the Asian option at a sequence of valuation dates.
#
# Before averaging begins, the option behaves similarly to a forward-starting
# path-dependent option. Once averaging has started, the accrued average
# becomes part of the state of the contract and influences the remaining
# option value.
#
# This example compares the different valuation methods as the valuation date
# moves towards expiry.
# ============================================================================

print("\n" + "=" * 78)
print("11. TIME EVOLUTION")
print("=" * 78)

start_averaging_dt_time = Date(1, 1, 2015)
expiry_dt_time = Date(1, 1, 2016)

stock_price_time = 100.0
volatility_time = 0.20
interest_rate_time = 0.30
dividend_yield_time = 0.10

num_obs_per_year_time = 100
strike_price_time = 100.0

accrued_average_time = stock_price_time * 0.90

model_time = BlackScholes(
    volatility_time,
)

asian_option_time = EquityAsianOption(
    start_averaging_dt_time,
    expiry_dt_time,
    strike_price_time,
    OptionTypes.EUROPEAN_CALL,
    num_obs_per_year_time,
)

value_dts = [
    Date(1, 4, 2014),
    Date(1, 6, 2014),
    Date(1, 8, 2014),
    Date(1, 2, 2015),
    Date(1, 4, 2015),
    Date(1, 6, 2015),
    Date(1, 8, 2015),
]

time_turnbull_values = []
time_curran_values = []
time_geometric_values = []
time_mc_values = []

num_paths_time = 10000

print(
    f"{'DATE':>15s}"
    f"{'GEOMETRIC':>18s}"
    f"{'TURNBULL':>18s}"
    f"{'CURRAN':>18s}"
    f"{'MONTE CARLO':>18s}"
)

print("-" * 87)

for time_value_dt in value_dts:

    time_discount_curve = FlatDiscountCurve(
        time_value_dt,
        interest_rate_time,
    )

    time_dividend_curve = FlatDiscountCurve(
        time_value_dt,
        dividend_yield_time,
    )

    geometric_value = asian_option_time.value(
        time_value_dt,
        stock_price_time,
        time_discount_curve,
        time_dividend_curve,
        model_time,
        AsianOptionValuationTypes.GEOMETRIC,
        accrued_average_time,
    )

    turnbull_value = asian_option_time.value(
        time_value_dt,
        stock_price_time,
        time_discount_curve,
        time_dividend_curve,
        model_time,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_average_time,
    )

    curran_value = asian_option_time.value(
        time_value_dt,
        stock_price_time,
        time_discount_curve,
        time_dividend_curve,
        model_time,
        AsianOptionValuationTypes.CURRAN,
        accrued_average_time,
    )

    mc_value = asian_option_time.value_mc(
        time_value_dt,
        stock_price_time,
        time_discount_curve,
        time_dividend_curve,
        model_time,
        num_paths_time,
        seed,
        accrued_average_time,
    )

    time_geometric_values.append(
        geometric_value,
    )

    time_turnbull_values.append(
        turnbull_value,
    )

    time_curran_values.append(
        curran_value,
    )

    time_mc_values.append(
        mc_value,
    )

    print(
        f"{str(time_value_dt):>15s}"
        f"{geometric_value:18.8f}"
        f"{turnbull_value:18.8f}"
        f"{curran_value:18.8f}"
        f"{mc_value:18.8f}"
    )


# ============================================================================
# 12. PLOT TIME EVOLUTION
# ============================================================================
# What this section demonstrates:
# Plots the value produced by each valuation method as the valuation date moves
# towards expiry.
#
# Differences between the methods can become particularly interesting after
# the averaging period has started because the realised average affects the
# remaining payoff distribution.
# ============================================================================

print("\n" + "=" * 78)
print("12. PLOT TIME EVOLUTION")
print("=" * 78)

# FinancePy Date objects are not Python datetime objects, so use integer
# positions on the x-axis and display the FinancePy dates as tick labels.
plot_dates = np.arange(len(value_dts))
plot_date_labels = [str(dt) for dt in value_dts]

plt.figure(figsize=(9, 6))

plt.plot(
    plot_dates,
    time_geometric_values,
    marker="o",
    label="Geometric",
)

plt.plot(
    plot_dates,
    time_turnbull_values,
    marker="o",
    label="Turnbull-Wakeman",
)

plt.plot(
    plot_dates,
    time_curran_values,
    marker="o",
    label="Curran",
)

plt.plot(
    plot_dates,
    time_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Valuation Date")
plt.ylabel("Asian Option Value")
plt.title("Asian Option Value through Time")

plt.xticks(
    plot_dates,
    plot_date_labels,
    rotation=45,
)

plt.grid(True)
plt.legend()
plt.show()

# =============================================================================
# MONTE CARLO IMPLEMENTATION CONSISTENCY
# =============================================================================
# The standard and fast Monte Carlo implementations consume random numbers
# in a different order. Therefore, for a given seed they need not produce
# the same value.
#
# This test runs both implementations over many independent seeds and checks
# whether the mean difference between them is statistically consistent with
# zero.

print("\n" + "=" * 78)
print("MONTE CARLO IMPLEMENTATION CONSISTENCY")
print("=" * 78)

########################################################################################
# Market / contract
value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
start_averaging_dt = value_dt

stock_price = 100.0
strike_price = 100.0

interest_rate = 0.05
dividend_yield = 0.01
volatility = 0.30

# Observation frequency
num_obs_per_year = 252

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)
model = BlackScholes(volatility)

option = EquityAsianOption(
    start_averaging_dt,
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    num_obs_per_year,
)

num_paths = 10000
num_replications = 100
seed_start = 1000

values_mc = []
values_mc_fast = []

for seed in range(seed_start, seed_start + num_replications):

    value_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths=num_paths,
        seed=seed,
        accrued_average=None,
    )

    value_mc_fast = option.value_mc_fast(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths=num_paths,
        seed=seed,
        accrued_average=None,
    )

    values_mc.append(value_mc)
    values_mc_fast.append(value_mc_fast)

values_mc = np.asarray(values_mc)
values_mc_fast = np.asarray(values_mc_fast)

diffs = values_mc - values_mc_fast

mean_mc = np.mean(values_mc)
mean_mc_fast = np.mean(values_mc_fast)

std_mc = np.std(values_mc, ddof=1)
std_mc_fast = np.std(values_mc_fast, ddof=1)

se_mc = std_mc / np.sqrt(num_replications)
se_mc_fast = std_mc_fast / np.sqrt(num_replications)

mean_diff = np.mean(diffs)
std_diff = np.std(diffs, ddof=1)
se_diff = std_diff / np.sqrt(num_replications)

if se_diff > 0.0:
    mean_over_se = mean_diff / se_diff
else:
    mean_over_se = 0.0

print(f"Paths             : {num_paths}")
print(f"Replications      : {num_replications}")
print()

print(
    f"{'METHOD':<18s}"
    f"{'MEAN':>14s}"
    f"{'STD':>14s}"
    f"{'MEAN SE':>14s}"
)

print("-" * 60)

print(
    f"{'MC':<18s}"
    f"{mean_mc:14.8f}"
    f"{std_mc:14.8f}"
    f"{se_mc:14.8f}"
)

print(
    f"{'MC Fast':<18s}"
    f"{mean_mc_fast:14.8f}"
    f"{std_mc_fast:14.8f}"
    f"{se_mc_fast:14.8f}"
)

print()
print(f"Mean difference   : {mean_diff:.10f}")
print(f"Std difference    : {std_diff:.10f}")
print(f"SE difference     : {se_diff:.10f}")
print(f"Mean / SE         : {mean_over_se:.4f}")


# =============================================================================
# 13. MONTE CARLO IMPLEMENTATION CONSISTENCY
# =============================================================================
# Compare the standard, fast, and control-variate Monte Carlo implementations
# over many independent seeds.
#
# Standard MC and Fast MC consume random numbers in different orders, so they
# are not expected to give identical values for an individual seed. Their
# replicated means, however, should be statistically consistent.
#
# The control-variate estimator should have substantially lower variance and
# therefore provides a more precise Monte Carlo benchmark.

print("\n" + "=" * 78)
print("13. MONTE CARLO IMPLEMENTATION CONSISTENCY")
print("=" * 78)

num_paths = 10000
num_replications = 100
seed_start = 1000

values_mc = []
values_mc_fast = []
values_mc_cv = []

for seed in range(seed_start, seed_start + num_replications):

    v_mc = option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths=num_paths,
        seed=seed,
        accrued_average=None,
    )

    v_mc_fast = option.value_mc_fast(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths=num_paths,
        seed=seed,
        accrued_average=None,
    )

    v_mc_cv = option.value_mc_fast_cv(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths=num_paths,
        seed=seed,
        accrued_average=None,
    )

    values_mc.append(v_mc)
    values_mc_fast.append(v_mc_fast)
    values_mc_cv.append(v_mc_cv)

values_mc = np.asarray(values_mc)
values_mc_fast = np.asarray(values_mc_fast)
values_mc_cv = np.asarray(values_mc_cv)


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

mean_mc = np.mean(values_mc)
mean_mc_fast = np.mean(values_mc_fast)
mean_mc_cv = np.mean(values_mc_cv)

std_mc = np.std(values_mc, ddof=1)
std_mc_fast = np.std(values_mc_fast, ddof=1)
std_mc_cv = np.std(values_mc_cv, ddof=1)

se_mc = std_mc / np.sqrt(num_replications)
se_mc_fast = std_mc_fast / np.sqrt(num_replications)
se_mc_cv = std_mc_cv / np.sqrt(num_replications)


# -----------------------------------------------------------------------------
# Standard MC versus Fast MC
# -----------------------------------------------------------------------------

diffs = values_mc - values_mc_fast

mean_diff = np.mean(diffs)
std_diff = np.std(diffs, ddof=1)
se_diff = std_diff / np.sqrt(num_replications)

if se_diff > 0.0:
    mean_over_se = mean_diff / se_diff
else:
    mean_over_se = 0.0


print(f"Paths             : {num_paths}")
print(f"Replications      : {num_replications}")
print()

print(
    f"{'METHOD':<18s}"
    f"{'MEAN':>14s}"
    f"{'STD':>14s}"
    f"{'MEAN SE':>14s}"
)

print("-" * 60)

print(
    f"{'MC':<18s}"
    f"{mean_mc:14.8f}"
    f"{std_mc:14.8f}"
    f"{se_mc:14.8f}"
)

print(
    f"{'MC Fast':<18s}"
    f"{mean_mc_fast:14.8f}"
    f"{std_mc_fast:14.8f}"
    f"{se_mc_fast:14.8f}"
)

print(
    f"{'MC Fast CV':<18s}"
    f"{mean_mc_cv:14.8f}"
    f"{std_mc_cv:14.8f}"
    f"{se_mc_cv:14.8f}"
)

print()
print("MC versus MC Fast")
print("-" * 40)
print(f"Mean difference   : {mean_diff:.10f}")
print(f"Std difference    : {std_diff:.10f}")
print(f"SE difference     : {se_diff:.10f}")
print(f"Mean / SE         : {mean_over_se:.4f}")


# -----------------------------------------------------------------------------
# Variance reduction
# -----------------------------------------------------------------------------

if std_mc_cv > 0.0:
    variance_reduction = (std_mc_fast / std_mc_cv) ** 2
else:
    variance_reduction = np.inf

print()
print("Control variate")
print("-" * 40)
print(f"Fast MC std       : {std_mc_fast:.10f}")
print(f"Fast MC CV std    : {std_mc_cv:.10f}")
print(f"Variance reduction: {variance_reduction:.2f}x")


# =============================================================================
# 14. ANALYTIC APPROXIMATIONS VERSUS MONTE CARLO
# =============================================================================
# Use the control-variate Monte Carlo estimator as the numerical benchmark.
# Curran and Turnbull-Wakeman are approximations, so differences from the
# Monte Carlo benchmark are expected and are not necessarily errors.

print("\n" + "=" * 78)
print("14. ANALYTIC APPROXIMATIONS VERSUS MONTE CARLO")
print("=" * 78)

value_curran = option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.CURRAN,
    accrued_average=None,
)

value_tw = option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    accrued_average=None,
)

# Use control-variate MC as the numerical benchmark.
mc_reference = mean_mc_cv
mc_reference_se = se_mc_cv

print(
    f"{'METHOD':<24s}"
    f"{'VALUE':>14s}"
    f"{'DIFF FROM CV MC':>18s}"
)

print("-" * 56)

print(
    f"{'MC':<24s}"
    f"{mean_mc:14.8f}"
    f"{mean_mc - mc_reference:18.8f}"
)

print(
    f"{'MC Fast':<24s}"
    f"{mean_mc_fast:14.8f}"
    f"{mean_mc_fast - mc_reference:18.8f}"
)

print(
    f"{'MC Fast CV':<24s}"
    f"{mean_mc_cv:14.8f}"
    f"{0.0:18.8f}"
)

print(
    f"{'Curran':<24s}"
    f"{value_curran:14.8f}"
    f"{value_curran - mc_reference:18.8f}"
)

print(
    f"{'Turnbull-Wakeman':<24s}"
    f"{value_tw:14.8f}"
    f"{value_tw - mc_reference:18.8f}"
)

print()
print(f"CV MC mean standard error : {mc_reference_se:.10f}")

# =============================================================================
# 15. OBSERVATION-FREQUENCY CONVERGENCE
# =============================================================================
# Compare the analytic approximations with control-variate Monte Carlo as
# the observation frequency increases.
#
# The control-variate Monte Carlo estimator is used as the numerical
# benchmark because it has substantially lower variance than standard MC.

print("\n" + "=" * 78)
print("15. OBSERVATION-FREQUENCY CONVERGENCE")
print("=" * 78)

num_obs_per_year_list = [
    12,
    26,
    52,
    100,
    252,
    500,
    1000,
    2000,
]


num_paths = 10000
num_replications = 50
seed_start = 2000

results = []

print(f"Paths             : {num_paths}")
print(f"Replications      : {num_replications}")
print()

print(
    f"{'OBS/YEAR':>10s}"
    f"{'CV MC':>14s}"
    f"{'CV SE':>12s}"
    f"{'CURRAN':>14s}"
    f"{'CURRAN-CV':>14s}"
    f"{'TW':>14s}"
    f"{'TW-CV':>14s}"
)

print("-" * 92)

for num_obs_per_year in num_obs_per_year_list:

    option = EquityAsianOption(
        start_averaging_dt,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
        num_obs_per_year,
    )

    # -------------------------------------------------------------------------
    # Replicated control-variate Monte Carlo
    # -------------------------------------------------------------------------

    cv_values = []

    for seed in range(seed_start, seed_start + num_replications):

        v_cv = option.value_mc_fast_cv(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths=num_paths,
            seed=seed,
            accrued_average=None,
        )

        cv_values.append(v_cv)

    cv_values = np.asarray(cv_values)

    mean_cv = np.mean(cv_values)
    std_cv = np.std(cv_values, ddof=1)
    se_cv = std_cv / np.sqrt(num_replications)

    # -------------------------------------------------------------------------
    # Curran
    # -------------------------------------------------------------------------

    v_curran = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.CURRAN,
        accrued_average=None,
    )

    # -------------------------------------------------------------------------
    # Turnbull-Wakeman
    # -------------------------------------------------------------------------

    v_tw = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_average=None,
    )

    curran_error = v_curran - mean_cv
    tw_error = v_tw - mean_cv

    results.append(
        (
            num_obs_per_year,
            mean_cv,
            se_cv,
            v_curran,
            curran_error,
            v_tw,
            tw_error,
        )
    )

    print(
        f"{num_obs_per_year:10d}"
        f"{mean_cv:14.8f}"
        f"{se_cv:12.8f}"
        f"{v_curran:14.8f}"
        f"{curran_error:14.8f}"
        f"{v_tw:14.8f}"
        f"{tw_error:14.8f}"
    )


# =============================================================================
# PLOT OBSERVATION-FREQUENCY CONVERGENCE
# =============================================================================

obs = np.array([x[0] for x in results])

cv_values = np.array([x[1] for x in results])
cv_se = np.array([x[2] for x in results])

curran_values = np.array([x[3] for x in results])
tw_values = np.array([x[5] for x in results])

plt.figure(figsize=(10, 6))

plt.errorbar(
    obs,
    cv_values,
    yerr=1.96 * cv_se,
    marker="o",
    capsize=3,
    label="MC Fast CV (95% CI)",
)

plt.plot(
    obs,
    curran_values,
    marker="o",
    label="Curran",
)

plt.plot(
    obs,
    tw_values,
    marker="o",
    label="Turnbull-Wakeman",
)

plt.xscale("log")

plt.xlabel("Observations per year")
plt.ylabel("Option value")
plt.title("Asian Option Value versus Observation Frequency")
plt.grid(True)
plt.legend()
plt.show()


print(
    f"{'OBS/YEAR':>10s}"
    f"{'MC FAST':>14s}"
    f"{'MC CV':>14s}"
    f"{'CV-FAST':>14s}"
    f"{'CURRAN':>14s}"
    f"{'TW':>14s}"
)

print("-" * 80)

for num_obs_per_year in [
    52,
    100,
    252,
    500,
    1000,
    2000,
]:

    option = EquityAsianOption(
        start_averaging_dt,
        expiry_dt,
        strike_price,
        OptionTypes.EUROPEAN_CALL,
        num_obs_per_year,
    )

    fast_values = []
    cv_values = []

    for seed in range(2000, 2050):

        v_fast = option.value_mc_fast(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths=10000,
            seed=seed,
            accrued_average=None,
        )

        v_cv = option.value_mc_fast_cv(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
            num_paths=10000,
            seed=seed,
            accrued_average=None,
        )

        fast_values.append(v_fast)
        cv_values.append(v_cv)

    mean_fast = np.mean(fast_values)
    mean_cv = np.mean(cv_values)

    v_curran = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.CURRAN,
        accrued_average=None,
    )

    v_tw = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        AsianOptionValuationTypes.TURNBULL_WAKEMAN,
        accrued_average=None,
    )

    print(
        f"{num_obs_per_year:10d}"
        f"{mean_fast:14.8f}"
        f"{mean_cv:14.8f}"
        f"{mean_cv - mean_fast:14.8f}"
        f"{v_curran:14.8f}"
        f"{v_tw:14.8f}"
    )
