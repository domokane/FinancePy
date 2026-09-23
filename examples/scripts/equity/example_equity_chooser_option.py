# ============================================================================
# FINANCEPY EXAMPLES - EquityChooserOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

from financepy.products.equity.equity_chooser_option import EquityChooserOption
from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


# ============================================================================
# 1. EQUITY CHOOSER OPTION HAUG
# ============================================================================
# What this section demonstrates:
# Checks the FinancePy analytical chooser-option value against the published
# Haug benchmark.
#
# It also compares the analytical value with Monte Carlo. The Monte Carlo
# result will not agree exactly because of simulation noise.

print("\n" + "=" * 78)
print("1. EQUITY CHOOSER OPTION HAUG")
print("=" * 78)

# Following example in Haug Page 130.

value_dt = Date(1, 1, 2015)
choose_dt = Date(2, 4, 2015)
call_expiry_dt = Date(1, 7, 2015)
put_expiry_dt = Date(2, 8, 2015)

call_strike = 55.0
put_strike = 48.0

stock_price = 50.0
volatility = 0.35
interest_rate = 0.10
dividend_yield = 0.05

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

chooser_option = EquityChooserOption(
    choose_dt,
    call_expiry_dt,
    put_expiry_dt,
    call_strike,
    put_strike,
)

v = chooser_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

num_paths = 100000

v_mc = chooser_option.value_mc(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
)

v_haug = 6.0508

print("METHOD                 VALUE          DIFF VS HAUG")
print("-" * 55)
print(f"FINANCEPY ANALYTIC {v:14.8f} {v - v_haug:18.8f}")
print(f"HAUG                {v_haug:14.8f} {0.0:18.8f}")
print(f"MONTE CARLO         {v_mc:14.8f} {v_mc - v_haug:18.8f}")


# ============================================================================
# 2. EQUITY CHOOSER OPTION MATLAB
# ============================================================================
# What this section demonstrates:
# Checks the FinancePy analytical chooser-option value against the MATLAB
# example.
#
# This provides an independent benchmark using equal call and put strikes
# and equal option expiry dates.

print("\n" + "=" * 78)
print("2. EQUITY CHOOSER OPTION MATLAB")
print("=" * 78)

# Reference:
# https://fr.mathworks.com/help/fininst/chooserbybls.html

value_dt = Date(1, 6, 2007)
choose_dt = Date(31, 8, 2007)
call_expiry_dt = Date(2, 12, 2007)
put_expiry_dt = Date(2, 12, 2007)

call_strike = 60.0
put_strike = 60.0

stock_price = 50.0
volatility = 0.20
interest_rate = 0.10
dividend_yield = 0.05

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

chooser_option = EquityChooserOption(
    choose_dt,
    call_expiry_dt,
    put_expiry_dt,
    call_strike,
    put_strike,
)

v = chooser_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

v_mc = chooser_option.value_mc(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
)

v_matlab = 8.9308

print("METHOD                 VALUE        DIFF VS MATLAB")
print("-" * 55)
print(f"FINANCEPY ANALYTIC {v:14.8f} {v - v_matlab:18.8f}")
print(f"MATLAB              {v_matlab:14.8f} {0.0:18.8f}")
print(f"MONTE CARLO         {v_mc:14.8f} {v_mc - v_matlab:18.8f}")


# ============================================================================
# 3. EQUITY CHOOSER OPTION DERIVICOM
# ============================================================================
# What this section demonstrates:
# Checks the FinancePy analytical chooser-option value against the Derivicom
# benchmark.
#
# Unlike the MATLAB example, the call and put have different strikes and
# different expiry dates.

print("\n" + "=" * 78)
print("3. EQUITY CHOOSER OPTION DERIVICOM")
print("=" * 78)

# Reference:
# http://derivicom.com/support/finoptionsxl/index.html?complex_chooser.htm

value_dt = Date(1, 1, 2007)
choose_dt = Date(1, 2, 2007)
call_expiry_dt = Date(1, 4, 2007)
put_expiry_dt = Date(1, 5, 2007)

call_strike = 40.0
put_strike = 35.0

stock_price = 38.0
volatility = 0.20
interest_rate = 0.08
dividend_yield = 0.0625

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

chooser_option = EquityChooserOption(
    choose_dt,
    call_expiry_dt,
    put_expiry_dt,
    call_strike,
    put_strike,
)

v = chooser_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

v_mc = chooser_option.value_mc(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
)

v_derivicom = 1.0989

print("METHOD                 VALUE     DIFF VS DERIVICOM")
print("-" * 55)
print(f"FINANCEPY ANALYTIC {v:14.8f} {v - v_derivicom:18.8f}")
print(f"DERIVICOM           {v_derivicom:14.8f} {0.0:18.8f}")
print(f"MONTE CARLO         {v_mc:14.8f} {v_mc - v_derivicom:18.8f}")


# ============================================================================
# 4. CHOOSER OPTION VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how the chooser-option value changes with the current stock price.
#
# A chooser option has value on both sides of the strike because the holder
# will later choose whether the contract becomes a call or a put.
#
# This produces the characteristic U-shaped chooser-option value profile.

print("\n" + "=" * 78)
print("4. CHOOSER OPTION VALUE VERSUS STOCK PRICE")
print("=" * 78)

value_dt = Date(1, 1, 2027)
choose_dt = Date(1, 6, 2027)

call_expiry_dt = Date(1, 1, 2028)
put_expiry_dt = Date(1, 1, 2028)

call_strike = 100.0
put_strike = 100.0

volatility = 0.20
interest_rate = 0.04
dividend_yield = 0.02

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

chooser_option = EquityChooserOption(
    choose_dt,
    call_expiry_dt,
    put_expiry_dt,
    call_strike,
    put_strike,
)

stock_prices = np.linspace(
    50.0,
    150.0,
    101,
)

chooser_values = chooser_option.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curve,
    model,
)

print("STOCK PRICE       CHOOSER VALUE")
print("-" * 35)

for stock_price, chooser_value in zip(
    stock_prices[::10],
    chooser_values[::10],
):
    print(
        f"{stock_price:11.2f}"
        f"{chooser_value:20.8f}"
    )

plt.figure()

plt.plot(
    stock_prices,
    chooser_values,
    label="Chooser Option",
)

plt.axvline(
    call_strike,
    linestyle="--",
    label="Strike",
)

plt.xlabel("Stock Price")
plt.ylabel("Chooser Option Value")
plt.title("Chooser Option Value versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. VALUE OF THE CHOICE FEATURE
# ============================================================================
# What this section demonstrates:
# Compares the chooser option with the corresponding European call and put.
#
# Buying a chooser option is not the same as choosing the more valuable
# vanilla option today. The holder waits until the future choice date before
# deciding whether to receive the call or the put.
#
# The difference therefore measures the value associated with retaining that
# future choice.

print("\n" + "=" * 78)
print("5. VALUE OF THE CHOICE FEATURE")
print("=" * 78)

call_option = EquityVanillaOption(
    call_expiry_dt,
    call_strike,
    OptionTypes.EUROPEAN_CALL,
)

put_option = EquityVanillaOption(
    put_expiry_dt,
    put_strike,
    OptionTypes.EUROPEAN_PUT,
)

call_values = np.zeros(len(stock_prices))
put_values = np.zeros(len(stock_prices))

for i, stock_price in enumerate(stock_prices):

    call_values[i] = call_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    put_values[i] = put_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

best_today_values = np.maximum(
    call_values,
    put_values,
)

choice_feature_values = (
    chooser_values
    - best_today_values
)

print(
    "STOCK"
    "        CALL"
    "         PUT"
    "     BEST TODAY"
    "       CHOOSER"
    "    CHOICE VALUE"
)

print("-" * 82)

for i in range(
    0,
    len(stock_prices),
    10,
):
    print(
        f"{stock_prices[i]:7.2f}"
        f"{call_values[i]:13.6f}"
        f"{put_values[i]:13.6f}"
        f"{best_today_values[i]:15.6f}"
        f"{chooser_values[i]:14.6f}"
        f"{choice_feature_values[i]:16.6f}"
    )

plt.figure()

plt.plot(
    stock_prices,
    chooser_values,
    label="Chooser",
)

plt.plot(
    stock_prices,
    call_values,
    label="European Call",
)

plt.plot(
    stock_prices,
    put_values,
    label="European Put",
)

plt.plot(
    stock_prices,
    best_today_values,
    linestyle="--",
    label="Best Vanilla Today",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("Chooser Option versus Vanilla Options")
plt.grid(True)
plt.legend()
plt.show()


plt.figure()

plt.plot(
    stock_prices,
    choice_feature_values,
    label="Value of Future Choice",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Stock Price")
plt.ylabel("Chooser - Best Vanilla Today")
plt.title("Value of Retaining the Future Choice")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. CHOOSER VALUE VERSUS CHOICE DATE
# ============================================================================
# What this section demonstrates:
# Shows how the chooser-option value changes when the date on which the holder
# must choose between the call and put is changed.
#
# A later choice date generally gives the holder more information before the
# decision has to be made, so this experiment illustrates the economic value
# of decision flexibility.

print("\n" + "=" * 78)
print("6. CHOOSER VALUE VERSUS CHOICE DATE")
print("=" * 78)

value_dt = Date(1, 1, 2027)

call_expiry_dt = Date(1, 1, 2028)
put_expiry_dt = Date(1, 1, 2028)

call_strike = 100.0
put_strike = 100.0

stock_price = 100.0
volatility = 0.20
interest_rate = 0.04
dividend_yield = 0.02

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

choice_dates = [
    Date(1, 2, 2027),
    Date(1, 3, 2027),
    Date(1, 4, 2027),
    Date(1, 5, 2027),
    Date(1, 6, 2027),
    Date(1, 7, 2027),
    Date(1, 8, 2027),
    Date(1, 9, 2027),
    Date(1, 10, 2027),
    Date(1, 11, 2027),
    Date(1, 12, 2027),
]

choice_times = []
choice_values = []

print("CHOICE DATE          TIME        VALUE")
print("-" * 45)

for choose_dt in choice_dates:

    chooser_option = EquityChooserOption(
        choose_dt,
        call_expiry_dt,
        put_expiry_dt,
        call_strike,
        put_strike,
    )

    chooser_value = chooser_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    choice_time = (
        choose_dt - value_dt
    ) / 365.0

    choice_times.append(
        choice_time,
    )

    choice_values.append(
        chooser_value,
    )

    print(
        f"{str(choose_dt):15s}"
        f"{choice_time:12.6f}"
        f"{chooser_value:14.8f}"
    )

plt.figure()

plt.plot(
    choice_times,
    choice_values,
    marker="o",
)

plt.xlabel("Time to Choice Date (Years)")
plt.ylabel("Chooser Option Value")
plt.title("Chooser Option Value versus Choice Date")
plt.grid(True)
plt.show()


# ============================================================================
# 7. MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Compares Monte Carlo chooser-option values with the analytical FinancePy
# result as the number of simulation paths increases.
#
# Monte Carlo estimates contain sampling noise. Increasing the number of
# paths should reduce that noise, although convergence will not be perfectly
# monotonic for a single sequence of simulations.

print("\n" + "=" * 78)
print("7. MONTE CARLO CONVERGENCE")
print("=" * 78)

value_dt = Date(1, 1, 2015)
choose_dt = Date(2, 4, 2015)

call_expiry_dt = Date(1, 7, 2015)
put_expiry_dt = Date(2, 8, 2015)

call_strike = 55.0
put_strike = 48.0

stock_price = 50.0
volatility = 0.35
interest_rate = 0.10
dividend_yield = 0.05

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

chooser_option = EquityChooserOption(
    choose_dt,
    call_expiry_dt,
    put_expiry_dt,
    call_strike,
    put_strike,
)

analytic_value = chooser_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

num_paths_list = np.array(
    [
        1000,
        2000,
        5000,
        10000,
        20000,
        50000,
        100000,
    ]
)

mc_values = []
mc_errors = []

print(
    "NUM PATHS"
    "        ANALYTIC"
    "              MC"
    "           ERROR"
)

print("-" * 60)

for num_paths in num_paths_list:

    mc_value = chooser_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        int(num_paths),
    )

    mc_error = (
        mc_value
        - analytic_value
    )

    mc_values.append(
        mc_value,
    )

    mc_errors.append(
        mc_error,
    )

    print(
        f"{num_paths:9d}"
        f"{analytic_value:16.8f}"
        f"{mc_value:16.8f}"
        f"{mc_error:16.8f}"
    )

mc_values = np.array(
    mc_values,
)

mc_errors = np.array(
    mc_errors,
)

plt.figure()

plt.plot(
    num_paths_list,
    mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.axhline(
    analytic_value,
    linestyle="--",
    label="Analytic",
)

plt.xscale("log")

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Chooser Option Value")
plt.title("Chooser Option Monte Carlo Convergence")
plt.grid(True)
plt.legend()
plt.show()


plt.figure()

plt.plot(
    num_paths_list,
    mc_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xscale("log")

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Monte Carlo - Analytic")
plt.title("Chooser Option Monte Carlo Error")
plt.grid(True)
plt.show()
