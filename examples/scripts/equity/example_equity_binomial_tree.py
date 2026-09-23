# ============================================================================
# FINANCEPY EXAMPLES - EquityBinomialTree
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

from financepy.models.black_scholes import BlackScholes

from financepy.products.equity.equity_binomial_tree import EquityBinomialTree
from financepy.products.equity.equity_binomial_tree import EquityTreeExerciseTypes
from financepy.products.equity.equity_binomial_tree import EquityTreePayoffTypes
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption

LINE = "=" * 78


# ============================================================================
# SUPPORTING FUNCTION
# ============================================================================


def tree_value(
    tree,
    stock_price,
    discount_curve,
    dividend_curve,
    volatility,
    num_steps,
    value_dt,
    expiry_dt,
    exercise_type,
    option_sign,
    strike_price,
):
    """Value a vanilla option with EquityBinomialTree."""

    payoff = EquityTreePayoffTypes.VANILLA_OPTION
    params = np.array([option_sign, strike_price])

    return tree.value(
        stock_price,
        discount_curve,
        dividend_curve,
        volatility,
        int(num_steps),
        value_dt,
        payoff,
        expiry_dt,
        payoff,
        exercise_type,
        params,
    )


# ============================================================================
# MARKET DATA
# ============================================================================

stock_price = 50.0
strike_price = 50.0
risk_free_rate = 0.06
dividend_yield = 0.04
volatility = 0.40

value_dt = Date(1, 1, 2016)
expiry_dt = Date(1, 1, 2017)

discount_curve = FlatDiscountCurve(value_dt, risk_free_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

model = BlackScholes(volatility)
tree = EquityBinomialTree()


# ============================================================================
# 1. EQUITY BINOMIAL TREE
# ============================================================================
# What this section demonstrates:
#
# Values European and American vanilla puts and calls using the binomial tree.
# European values are compared with Black-Scholes. American options can be
# exercised before expiry, so their exercise decision is handled at each node.
# ============================================================================

print("\n" + LINE)
print("1. EQUITY BINOMIAL TREE")
print(LINE)

num_steps_list = [100, 500, 1000]

for option_name, option_type, option_sign in [
    ("EUROPEAN PUT", OptionTypes.EUROPEAN_PUT, -1.0),
    ("EUROPEAN CALL", OptionTypes.EUROPEAN_CALL, 1.0),
]:
    print(f"\n================== {option_name} ==================")

    option = EquityVanillaOption(expiry_dt, strike_price, option_type)
    bs_value = option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
    bs_delta = option.delta(value_dt, stock_price, discount_curve, dividend_curve, model)
    bs_gamma = option.gamma(value_dt, stock_price, discount_curve, dividend_curve, model)
    bs_theta = option.theta(value_dt, stock_price, discount_curve, dividend_curve, model)

    print(f"{'BS VALUE':>15}{'BS DELTA':>15}{'BS GAMMA':>15}{'BS THETA':>15}")
    print(f"{bs_value:15.6f}{bs_delta:15.6f}{bs_gamma:15.6f}{bs_theta:15.6f}")
    print(f"{'STEPS':>10}{'VALUE':>15}{'DELTA':>15}{'GAMMA':>15}{'THETA':>15}{'TIME':>15}")
    print("-" * 85)

    for num_steps in num_steps_list:
        start = time.perf_counter()
        results = tree_value(
            tree, stock_price, discount_curve, dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN,
            option_sign, strike_price,
        )
        duration = time.perf_counter() - start
        print(
            f"{num_steps:10d}{results[0]:15.6f}{results[1]:15.6f}"
            f"{results[2]:15.6f}{results[3]:15.6f}{duration:15.6f}"
        )

for option_name, option_sign in [
    ("AMERICAN PUT", -1.0),
    ("AMERICAN CALL", 1.0),
]:
    print(f"\n================== {option_name} ==================")
    print(f"{'STEPS':>10}{'VALUE':>15}{'DELTA':>15}{'GAMMA':>15}{'THETA':>15}{'TIME':>15}")
    print("-" * 85)

    for num_steps in num_steps_list:
        start = time.perf_counter()
        results = tree_value(
            tree, stock_price, discount_curve, dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN,
            option_sign, strike_price,
        )
        duration = time.perf_counter() - start
        print(
            f"{num_steps:10d}{results[0]:15.6f}{results[1]:15.6f}"
            f"{results[2]:15.6f}{results[3]:15.6f}{duration:15.6f}"
        )


# ============================================================================
# 2. EUROPEAN OPTION TREE CONVERGENCE
# ============================================================================
# What this section demonstrates:
#
# European vanilla options provide an exact Black-Scholes benchmark for the
# tree. As the number of time steps increases, tree values should approach the
# corresponding Black-Scholes values. Convergence can oscillate.
# ============================================================================

print("\n" + LINE)
print("2. EUROPEAN OPTION TREE CONVERGENCE")
print(LINE)

num_steps_grid = np.array(
    [
        10,
        11,
        20,
        21,
        50,
        51,
        100,
        101,
        200,
        201,
        500,
        501,
        1000,
        1001,
    ]
)

put_option = EquityVanillaOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_PUT)
call_option = EquityVanillaOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL)

bs_put_value = put_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
bs_call_value = call_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)

put_tree_values = []
call_tree_values = []

for num_steps in num_steps_grid:
    put_tree_values.append(
        tree_value(
            tree, stock_price, discount_curve, dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN,
            -1.0, strike_price,
        )[0]
    )
    call_tree_values.append(
        tree_value(
            tree, stock_price, discount_curve, dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN,
            1.0, strike_price,
        )[0]
    )

put_tree_values = np.array(put_tree_values)
call_tree_values = np.array(call_tree_values)

print(f"{'STEPS':>10}{'PUT TREE':>15}{'PUT BS':>15}{'CALL TREE':>15}{'CALL BS':>15}")
print("-" * 70)
for i, num_steps in enumerate(num_steps_grid):
    print(
        f"{num_steps:10d}{put_tree_values[i]:15.6f}{bs_put_value:15.6f}"
        f"{call_tree_values[i]:15.6f}{bs_call_value:15.6f}"
    )


# ============================================================================
# 3. PLOT EUROPEAN TREE CONVERGENCE
# ============================================================================
# What this section demonstrates:
#
# Combines put and call convergence in one graph. The dashed horizontal lines
# are the Black-Scholes benchmarks.
# ============================================================================

print("\n" + LINE)
print("3. PLOT EUROPEAN TREE CONVERGENCE")
print(LINE)

plt.figure()
plt.plot(num_steps_grid, put_tree_values, marker="o", label="Put Tree")
plt.plot(num_steps_grid, call_tree_values, marker="o", label="Call Tree")
plt.axhline(bs_put_value, linestyle="--", label="Put Black-Scholes")
plt.axhline(bs_call_value, linestyle="--", label="Call Black-Scholes")
plt.xlabel("Number of Tree Steps")
plt.ylabel("Option Value")
plt.title("European Option Binomial Tree Convergence")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. TREE CONVERGENCE ERROR
# ============================================================================
# What this section demonstrates:
#
# Subtracting the Black-Scholes benchmark isolates the tree discretisation
# error. The errors should tend towards zero as the tree is refined.
# ============================================================================

print("\n" + LINE)
print("4. TREE CONVERGENCE ERROR")
print(LINE)

put_errors = put_tree_values - bs_put_value
call_errors = call_tree_values - bs_call_value

print(f"{'STEPS':>10}{'PUT ERROR':>20}{'CALL ERROR':>20}")
print("-" * 50)
for i, num_steps in enumerate(num_steps_grid):
    print(f"{num_steps:10d}{put_errors[i]:20.8f}{call_errors[i]:20.8f}")

plt.figure()
plt.plot(num_steps_grid, put_errors, marker="o", label="Put")
plt.plot(num_steps_grid, call_errors, marker="o", label="Call")
plt.axhline(0.0, linestyle="--")
plt.xlabel("Number of Tree Steps")
plt.ylabel("Tree Value - Black-Scholes Value")
plt.title("Binomial Tree Convergence Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. AMERICAN PUT EARLY-EXERCISE PREMIUM
# ============================================================================
# What this section demonstrates:
#
# An American option contains all European exercise opportunities plus the
# right to exercise early. Therefore American value should not be below the
# corresponding European value. Their difference is the early-exercise
# premium.
# ============================================================================

print("\n" + LINE)
print("5. AMERICAN PUT EARLY-EXERCISE PREMIUM")
print(LINE)

stock_grid = np.linspace(25.0, 75.0, 51)
num_steps = 1000

european_put_values = []
american_put_values = []

for stock in stock_grid:
    european_put_values.append(
        tree_value(
            tree, stock, discount_curve, dividend_curve, volatility, num_steps,
            value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN, -1.0,
            strike_price,
        )[0]
    )
    american_put_values.append(
        tree_value(
            tree, stock, discount_curve, dividend_curve, volatility, num_steps,
            value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN, -1.0,
            strike_price,
        )[0]
    )

european_put_values = np.array(european_put_values)
american_put_values = np.array(american_put_values)
put_exercise_premium = american_put_values - european_put_values

plt.figure()
plt.plot(stock_grid, european_put_values, label="European Put")
plt.plot(stock_grid, american_put_values, label="American Put")
plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("European and American Put Values")
plt.grid(True)
plt.legend()
plt.show()

plt.figure()
plt.plot(stock_grid, put_exercise_premium)
plt.axhline(0.0, linestyle="--")
plt.xlabel("Stock Price")
plt.ylabel("American - European Value")
plt.title("American Put Early-Exercise Premium")
plt.grid(True)
plt.show()

minimum_put_premium = np.min(put_exercise_premium)
tolerance = 1.0e-10

print(f"{'MINIMUM EARLY-EXERCISE PREMIUM':<40}{minimum_put_premium:15.10f}")
print("CHECK:", "PASSED" if minimum_put_premium >= -tolerance else "FAILED")


# ============================================================================
# 6. AMERICAN CALL EARLY-EXERCISE PREMIUM
# ============================================================================
# What this section demonstrates:
#
# With a positive dividend yield, early exercise of an American call can have
# value because exercising allows the holder to own the stock and receive its
# dividends.
# ============================================================================

print("\n" + LINE)
print("6. AMERICAN CALL EARLY-EXERCISE PREMIUM")
print(LINE)

european_call_values = []
american_call_values = []

for stock in stock_grid:
    european_call_values.append(
        tree_value(
            tree, stock, discount_curve, dividend_curve, volatility, num_steps,
            value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN, 1.0,
            strike_price,
        )[0]
    )
    american_call_values.append(
        tree_value(
            tree, stock, discount_curve, dividend_curve, volatility, num_steps,
            value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN, 1.0,
            strike_price,
        )[0]
    )

european_call_values = np.array(european_call_values)
american_call_values = np.array(american_call_values)
call_exercise_premium = american_call_values - european_call_values

plt.figure()
plt.plot(stock_grid, call_exercise_premium)
plt.axhline(0.0, linestyle="--")
plt.xlabel("Stock Price")
plt.ylabel("American - European Value")
plt.title("American Call Early-Exercise Premium")
plt.grid(True)
plt.show()

print(f"{'MAXIMUM EARLY-EXERCISE PREMIUM':<40}{np.max(call_exercise_premium):15.10f}")


# ============================================================================
# 7. ZERO-DIVIDEND AMERICAN CALL TEST
# ============================================================================
# What this section demonstrates:
#
# For a non-dividend-paying stock with a non-negative interest rate, early
# exercise of a vanilla American call should not add value. The American and
# European call values should therefore agree up to tree discretisation error.
# ============================================================================

print("\n" + LINE)
print("7. ZERO-DIVIDEND AMERICAN CALL TEST")
print(LINE)

zero_dividend_curve = FlatDiscountCurve(value_dt, 0.0)
zero_div_european_values = []
zero_div_american_values = []

for stock in stock_grid:
    zero_div_european_values.append(
        tree_value(
            tree, stock, discount_curve, zero_dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.EUROPEAN,
            1.0, strike_price,
        )[0]
    )
    zero_div_american_values.append(
        tree_value(
            tree, stock, discount_curve, zero_dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN,
            1.0, strike_price,
        )[0]
    )

zero_div_european_values = np.array(zero_div_european_values)
zero_div_american_values = np.array(zero_div_american_values)
zero_div_call_premium = zero_div_american_values - zero_div_european_values
max_abs_zero_div_premium = np.max(np.abs(zero_div_call_premium))

print(f"{'MAX |AMERICAN - EUROPEAN|':<40}{max_abs_zero_div_premium:15.10f}")

# The comparison uses the same tree discretisation for both contracts, so the
# difference should be very small. A modest tolerance avoids making the example
# dependent on tiny implementation/platform differences.
zero_div_tolerance = 1.0e-8
print("CHECK:", "PASSED" if max_abs_zero_div_premium <= zero_div_tolerance else "CHECK NUMERICAL TOLERANCE")

plt.figure()
plt.plot(stock_grid, zero_div_call_premium)
plt.axhline(0.0, linestyle="--")
plt.xlabel("Stock Price")
plt.ylabel("American - European Value")
plt.title("Zero-Dividend American Call Early-Exercise Premium")
plt.grid(True)
plt.show()


# ============================================================================
# 8. TREE COMPUTATION TIME
# ============================================================================
# What this section demonstrates:
#
# Increasing the number of tree steps improves time discretisation but also
# increases computational cost. Each size is timed repeatedly and the median
# is reported to reduce noise from a single timing observation.
# ============================================================================

print("\n" + LINE)
print("8. TREE COMPUTATION TIME")
print(LINE)

timing_steps = np.array([50, 100, 200, 500, 1000, 2000])
num_timing_runs = 5
timing_values = []

# Warm up once so first-call overhead is not included in the timing study.
tree_value(
    tree, stock_price, discount_curve, dividend_curve, volatility, 50,
    value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN, -1.0, strike_price,
)

for num_steps in timing_steps:
    run_times = []

    for _ in range(num_timing_runs):
        start = time.perf_counter()
        tree_value(
            tree, stock_price, discount_curve, dividend_curve, volatility,
            num_steps, value_dt, expiry_dt, EquityTreeExerciseTypes.AMERICAN,
            -1.0, strike_price,
        )
        run_times.append(time.perf_counter() - start)

    timing_values.append(np.median(run_times))

timing_values = np.array(timing_values)

print(f"{'STEPS':>12}{'MEDIAN TIME (SEC)':>22}")
print("-" * 34)
for num_steps, elapsed in zip(timing_steps, timing_values):
    print(f"{num_steps:12d}{elapsed:22.8f}")

plt.figure()
plt.plot(timing_steps, timing_values, marker="o")
plt.xlabel("Number of Tree Steps")
plt.ylabel("Median Execution Time (seconds)")
plt.title("Binomial Tree Computation Time")
plt.grid(True)
plt.show()
