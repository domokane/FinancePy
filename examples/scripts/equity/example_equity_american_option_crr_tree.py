# ============================================================================
# FINANCEPY EXAMPLES - EquityAmericanOption CRR Tree
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#


import time
import numpy as np
import matplotlib.pyplot as plt


from financepy.products.equity.equity_american_option import EquityAmericanOption
from financepy.utils.global_types import OptionTypes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

set_plot_style()


# ============================================================================
# 1. EQUITY AMERICAN OPTION
# ============================================================================
# What this section demonstrates:
# Values European and American calls and puts using the CRR tree.
# Compares option value and Greeks using the same market assumptions.
# The American option includes the right to exercise before expiry.

print("\n" + "=" * 78)
print("1. EQUITY AMERICAN OPTION")
print("=" * 78)

value_dt = Date(1, 1, 2016)
expiry_dt = Date(1, 1, 2017)
stock_price = 50.0
interest_rate = 0.06
dividend_yield = 0.04
volatility = 0.40
strike_price = 50.0

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

num_steps = 1000

model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    num_steps,
)

print("================== EUROPEAN PUT =======================")

put_option = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_PUT,
)

value = put_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

delta = put_option.delta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

gamma = put_option.gamma(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

theta = put_option.theta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

print("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
print(
    "EUROPEAN_PUT_TREE",
    value,
    delta,
    gamma,
    theta,
)

print("================== AMERICAN PUT =======================")

put_option = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.AMERICAN_PUT,
)

value = put_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

delta = put_option.delta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

gamma = put_option.gamma(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

theta = put_option.theta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

print("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
print(
    "AMERICAN_PUT_TREE",
    value,
    delta,
    gamma,
    theta,
)

print("================== EUROPEAN CALL =======================")

call_option = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
)

value = call_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

delta = call_option.delta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

gamma = call_option.gamma(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

theta = call_option.theta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

print("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
print(
    "EUROPEAN_CALL_TREE",
    value,
    delta,
    gamma,
    theta,
)

print("================== AMERICAN CALL =======================")

call_option = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.AMERICAN_CALL,
)

value = call_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

delta = call_option.delta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

gamma = call_option.gamma(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

theta = call_option.theta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

print("opt_type", "VALUE", "DELTA", "GAMMA", "THETA")
print(
    "AMERICAN_CALL_TREE",
    value,
    delta,
    gamma,
    theta,
)


# ============================================================================
# 2. CRR TREE CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Shows how European and American option values change as the number of CRR
# tree steps is increased.
# A sufficiently fine tree should approach a stable option value.

print("\n" + "=" * 78)
print("2. CRR TREE CONVERGENCE")
print("=" * 78)

european_put = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_PUT,
)

american_put = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.AMERICAN_PUT,
)

european_call = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
)

american_call = EquityAmericanOption(
    expiry_dt,
    strike_price,
    OptionTypes.AMERICAN_CALL,
)

num_steps_list = [
    10,
    20,
    50,
    100,
    200,
    500,
    1000,
    2000,
]

european_put_values = []
american_put_values = []
european_call_values = []
american_call_values = []

print(
    f"{'STEPS':>10}"
    f"{'EUR PUT':>15}"
    f"{'AMER PUT':>15}"
    f"{'EUR CALL':>15}"
    f"{'AMER CALL':>15}"
    f"{'TIME':>15}"
)

print("-" * 85)

for num_steps in num_steps_list:

    model = BlackScholes(
        volatility,
        BlackScholesTypes.CRR_TREE,
        num_steps,
    )

    start = time.time()

    european_put_value = european_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    american_put_value = american_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    european_call_value = european_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    american_call_value = american_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    end = time.time()

    duration = end - start

    european_put_values.append(
        european_put_value,
    )

    american_put_values.append(
        american_put_value,
    )

    european_call_values.append(
        european_call_value,
    )

    american_call_values.append(
        american_call_value,
    )

    print(
        f"{num_steps:10d}"
        f"{european_put_value:15.6f}"
        f"{american_put_value:15.6f}"
        f"{european_call_value:15.6f}"
        f"{american_call_value:15.6f}"
        f"{duration:15.6f}"
    )


# ============================================================================
# 3. CRR TREE CONVERGENCE ERROR
# ============================================================================
# What this section demonstrates:
# Measures numerical convergence of the CRR tree against a high-resolution
# tree benchmark.
#
# Plotting the raw option values hides the convergence behaviour because the
# values are already very close. Plotting the absolute pricing error makes the
# remaining numerical error visible.
# ============================================================================

print("\n" + "=" * 78)
print("3. CRR TREE CONVERGENCE ERROR")
print("=" * 78)

benchmark_steps = 10000

benchmark_model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    benchmark_steps,
)

benchmark_european_put = european_put.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

benchmark_american_put = american_put.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

benchmark_european_call = european_call.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

benchmark_american_call = american_call.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

european_put_errors = np.abs(
    np.asarray(european_put_values) - benchmark_european_put
)

american_put_errors = np.abs(
    np.asarray(american_put_values) - benchmark_american_put
)

european_call_errors = np.abs(
    np.asarray(european_call_values) - benchmark_european_call
)

american_call_errors = np.abs(
    np.asarray(american_call_values) - benchmark_american_call
)

print(
    f"{'STEPS':>10}"
    f"{'EUR PUT ERR':>18}"
    f"{'AMER PUT ERR':>18}"
    f"{'EUR CALL ERR':>18}"
    f"{'AMER CALL ERR':>18}"
)

print("-" * 82)

for i, num_steps in enumerate(num_steps_list):

    print(
        f"{num_steps:10d}"
        f"{european_put_errors[i]:18.10f}"
        f"{american_put_errors[i]:18.10f}"
        f"{european_call_errors[i]:18.10f}"
        f"{american_call_errors[i]:18.10f}"
    )

plt.figure(figsize=(9, 6))

plt.loglog(
    num_steps_list,
    european_put_errors,
    marker="o",
    label="European Put",
)

plt.loglog(
    num_steps_list,
    american_put_errors,
    marker="o",
    label="American Put",
)

plt.loglog(
    num_steps_list,
    european_call_errors,
    marker="o",
    label="European Call",
)

plt.loglog(
    num_steps_list,
    american_call_errors,
    marker="o",
    label="American Call",
)

plt.xlabel("Number of CRR Tree Steps")
plt.ylabel("Absolute Pricing Error")
plt.title("CRR Tree Convergence Error")

plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 4. AMERICAN EARLY EXERCISE PREMIUM
# ============================================================================
# What this section demonstrates:
# Compares American and European option values over a range of stock prices.
# The difference between the two values measures the value of the additional
# right to exercise the American option before expiry.

print("\n" + "=" * 78)
print("4. AMERICAN EARLY EXERCISE PREMIUM")
print("=" * 78)

num_steps = 1000

model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    num_steps,
)

stock_prices = np.linspace(
    20.0,
    80.0,
    31,
)

european_put_values = []
american_put_values = []

european_call_values = []
american_call_values = []

put_exercise_premiums = []
call_exercise_premiums = []

print(
    f"{'STOCK':>10}"
    f"{'EUR PUT':>15}"
    f"{'AMER PUT':>15}"
    f"{'PUT PREMIUM':>15}"
    f"{'EUR CALL':>15}"
    f"{'AMER CALL':>15}"
    f"{'CALL PREMIUM':>15}"
)

print("-" * 100)

for stock_price_i in stock_prices:

    european_put_value = european_put.value(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
        model,
    )

    american_put_value = american_put.value(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
        model,
    )

    european_call_value = european_call.value(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
        model,
    )

    american_call_value = american_call.value(
        value_dt,
        stock_price_i,
        discount_curve,
        dividend_curve,
        model,
    )

    put_exercise_premium = (
        american_put_value
        - european_put_value
    )

    call_exercise_premium = (
        american_call_value
        - european_call_value
    )

    european_put_values.append(
        european_put_value,
    )

    american_put_values.append(
        american_put_value,
    )

    european_call_values.append(
        european_call_value,
    )

    american_call_values.append(
        american_call_value,
    )

    put_exercise_premiums.append(
        put_exercise_premium,
    )

    call_exercise_premiums.append(
        call_exercise_premium,
    )

    print(
        f"{stock_price_i:10.2f}"
        f"{european_put_value:15.6f}"
        f"{american_put_value:15.6f}"
        f"{put_exercise_premium:15.6f}"
        f"{european_call_value:15.6f}"
        f"{american_call_value:15.6f}"
        f"{call_exercise_premium:15.6f}"
    )


# ============================================================================
# 5. EARLY EXERCISE PREMIUM PLOT
# ============================================================================
# What this section demonstrates:
# Plots American minus European option value against stock price.
# This isolates the value of the American early-exercise feature for calls
# and puts.

print("\n" + "=" * 78)
print("5. EARLY EXERCISE PREMIUM PLOT")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    put_exercise_premiums,
    label="Put Early Exercise Premium",
)

plt.plot(
    stock_prices,
    call_exercise_premiums,
    label="Call Early Exercise Premium",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Stock Price")
plt.ylabel("American Value - European Value")
plt.title("American Option Early Exercise Premium")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. EUROPEAN AND AMERICAN PUT VALUES
# ============================================================================
# What this section demonstrates:
# Compares the full European and American put values across stock prices.
# The American put is worth at least as much as the corresponding European
# put because it includes the additional right to exercise before expiry.

print("\n" + "=" * 78)
print("6. EUROPEAN AND AMERICAN PUT VALUES")
print("=" * 78)

put_intrinsic_values = np.maximum(
    strike_price - stock_prices,
    0.0,
)

plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    european_put_values,
    label="European Put",
)

plt.plot(
    stock_prices,
    american_put_values,
    label="American Put",
)

plt.plot(
    stock_prices,
    put_intrinsic_values,
    linestyle="--",
    label="Intrinsic Value",
)

plt.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("European and American Put Values")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. EUROPEAN AND AMERICAN CALL VALUES
# ============================================================================
# What this section demonstrates:
# Compares the full European and American call values across stock prices.
# The non-zero dividend yield is important because dividends can make early
# exercise of an American call economically valuable.

print("\n" + "=" * 78)
print("7. EUROPEAN AND AMERICAN CALL VALUES")
print("=" * 78)

call_intrinsic_values = np.maximum(
    stock_prices - strike_price,
    0.0,
)

plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    european_call_values,
    label="European Call",
)

plt.plot(
    stock_prices,
    american_call_values,
    label="American Call",
)

plt.plot(
    stock_prices,
    call_intrinsic_values,
    linestyle="--",
    label="Intrinsic Value",
)

plt.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("European and American Call Values")

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. DIVIDEND YIELD AND AMERICAN CALL VALUE
# ============================================================================
# What this section demonstrates:
# Shows how the value of the American early-exercise feature changes with the
# dividend yield.
# A call on a non-dividend-paying stock generally has no benefit from early
# exercise, while dividends can make early exercise valuable.

print("\n" + "=" * 78)
print("8. DIVIDEND YIELD AND AMERICAN CALL VALUE")
print("=" * 78)

dividend_yields = np.linspace(
    0.0,
    0.12,
    25,
)

european_call_dividend_values = []
american_call_dividend_values = []
call_dividend_premiums = []

print(
    f"{'DIV YIELD':>15}"
    f"{'EUR CALL':>15}"
    f"{'AMER CALL':>15}"
    f"{'EARLY EX PREM':>20}"
)

print("-" * 65)

for dividend_yield_i in dividend_yields:

    dividend_curve_i = FlatDiscountCurve(
        value_dt,
        dividend_yield_i,
    )

    european_call_value = european_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve_i,
        model,
    )

    american_call_value = american_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve_i,
        model,
    )

    exercise_premium = (
        american_call_value
        - european_call_value
    )

    european_call_dividend_values.append(
        european_call_value,
    )

    american_call_dividend_values.append(
        american_call_value,
    )

    call_dividend_premiums.append(
        exercise_premium,
    )

    print(
        f"{dividend_yield_i * 100.0:15.4f}"
        f"{european_call_value:15.6f}"
        f"{american_call_value:15.6f}"
        f"{exercise_premium:20.6f}"
    )


# ============================================================================
# 9. DIVIDEND YIELD AND EARLY EXERCISE PREMIUM PLOT
# ============================================================================
# What this section demonstrates:
# Plots the American call early-exercise premium against dividend yield.
# This illustrates why the dividend assumption matters when valuing an
# American call.

print("\n" + "=" * 78)
print("9. DIVIDEND YIELD AND EARLY EXERCISE PREMIUM PLOT")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    dividend_yields * 100.0,
    call_dividend_premiums,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Dividend Yield (%)")
plt.ylabel("American Call - European Call")
plt.title("Dividend Yield and American Call Early Exercise Premium")

plt.grid(True)
plt.show()
