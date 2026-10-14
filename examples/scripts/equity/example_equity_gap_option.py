# ============================================================================
# FINANCEPY EXAMPLES - EquityGapOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.utils.helpers import option_years

from financepy.models.black_scholes import BlackScholes

from financepy.products.equity.equity_gap_option import EquityGapOption

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


LINE = "=" * 78


# ============================================================================
# 1. EQUITY GAP OPTION SETUP
# ============================================================================
# What this section demonstrates:
#
# Sets up a European equity gap call.
#
# A gap option has two different price levels:
#
#       trigger price = determines whether the payoff is activated
#       payoff strike = determines the amount paid at expiry
#
# For a gap call:
#
#       payoff = (S(T) - K2) * I(S(T) > K1)
#
# where K1 is the trigger price and K2 is the payoff strike.
# ============================================================================

print("\n" + LINE)
print("1. EQUITY GAP OPTION SETUP")
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

trigger_price = 100.0
payoff_strike = 95.0

gap_option = EquityGapOption(
    expiry_dt,
    trigger_price,
    payoff_strike,
    OptionTypes.EUROPEAN_CALL,
)

print(gap_option)


# ============================================================================
# 2. GAP CALL VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
#
# Shows how changing the payoff strike changes the behaviour of a gap call
# while keeping the trigger price fixed.
#
# If K2 < K1, exercise always produces a positive payoff.
#
# If K2 = K1, the gap call reduces to a standard European call.
#
# If K2 > K1, the option can produce a negative payoff when the terminal
# stock price lies between the trigger price and the payoff strike.
# ============================================================================

print("\n" + LINE)
print("2. GAP CALL VALUE VERSUS STOCK PRICE")
print(LINE)

trigger_price = 100.0


payoff_strikes = [
    80.0,    # Positive digital component
    100.0,   # Ordinary European call
    120.0,   # Negative digital component
]

stock_grid = np.linspace(
    50.0,
    150.0,
    101,
)

plt.figure()

for payoff_strike in payoff_strikes:

    gap_call = EquityGapOption(
        expiry_dt,
        trigger_price,
        payoff_strike,
        OptionTypes.EUROPEAN_CALL,
    )

    option_values = []

    for stock_price in stock_grid:

        option_value = gap_call.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        option_values.append(
            option_value,
        )

    option_values = np.asarray(
        option_values,
    )

    plt.plot(
        stock_grid,
        option_values,
        label=f"Payoff Strike = {payoff_strike:.0f}",
    )

plt.axvline(
    trigger_price,
    linestyle="--",
    linewidth=1.0,
    label=f"Trigger Price = {trigger_price:.0f}",
)

plt.axhline(
    0.0,
    linestyle=":",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Gap Call Value versus Stock Price"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. GAP CALL DELTA VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
#
# Shows how changing the payoff strike changes the delta of a gap call while
# keeping the trigger price fixed.
#
# When the trigger price and payoff strike differ, the gap option has a
# discontinuous terminal payoff at the trigger price.
#
# This can produce substantially different hedge behaviour from a standard
# European call.
# ============================================================================

print("\n" + LINE)
print("3. GAP CALL DELTA VERSUS STOCK PRICE")
print(LINE)

plt.figure()

for payoff_strike in payoff_strikes:

    gap_call = EquityGapOption(
        expiry_dt,
        trigger_price,
        payoff_strike,
        OptionTypes.EUROPEAN_CALL,
    )

    delta_values = []

    for stock_price in stock_grid:

        delta = gap_call.delta(
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

    plt.plot(
        stock_grid,
        delta_values,
        label=f"Payoff Strike = {payoff_strike:.0f}",
    )

plt.axvline(
    trigger_price,
    linestyle="--",
    linewidth=1.0,
    label=f"Trigger Price = {trigger_price:.0f}",
)

plt.axhline(
    0.0,
    linestyle=":",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Delta")

plt.title(
    "Gap Call Delta versus Stock Price"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. GAP CALL VALUE NEAR EXPIRY
# ============================================================================
# What this section demonstrates:
#
# Shows how the value of a gap call changes as the time to expiry decreases.
#
# The trigger price is above the payoff strike:
#
#       trigger price K1 = 100
#       payoff strike K2 = 80
#
# The terminal payoff is therefore:
#
#       payoff = 0             if S(T) <= 100
#       payoff = S(T) - 80     if S(T) > 100
#
# Immediately above the trigger price, the terminal payoff is close to 20.
#
# As expiry approaches, the option value approaches this discontinuous
# terminal payoff.
# ============================================================================

print("\n" + LINE)
print("4. GAP CALL VALUE NEAR EXPIRY")
print(LINE)

value_dt = Date(1, 1, 2015)

trigger_price = 100.0
payoff_strike = 80.0

expiry_dates = [
    (value_dt.add_months(12), "1 Year"),
    (value_dt.add_months(6), "6 Months"),
    (value_dt.add_months(3), "3 Months"),
    (value_dt.add_months(1), "1 Month"),
    (value_dt.add_days(7), "1 Week"),
]

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

stock_grid = np.linspace(
    80.0,
    130.0,
    201,
)

plt.figure()

for expiry_dt_near, maturity_label in expiry_dates:

    gap_call = EquityGapOption(
        expiry_dt_near,
        trigger_price,
        payoff_strike,
        OptionTypes.EUROPEAN_CALL,
    )

    option_values = []

    for stock_price in stock_grid:

        option_value = gap_call.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )

        option_values.append(
            option_value,
        )

    option_values = np.asarray(
        option_values,
    )

    plt.plot(
        stock_grid,
        option_values,
        label=maturity_label,
    )

plt.axvline(
    trigger_price,
    linestyle="--",
    linewidth=1.0,
    label=f"Trigger Price = {trigger_price:.0f}",
)

plt.axvline(
    payoff_strike,
    linestyle=":",
    linewidth=1.0,
    label=f"Payoff Strike = {payoff_strike:.0f}",
)

plt.axhline(
    0.0,
    linestyle="--",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Gap Call Value Near Expiry"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. GAP CALL DELTA NEAR EXPIRY
# ============================================================================
# What this section demonstrates:
#
# Shows how the delta of a gap call changes as time to expiry decreases.
#
# Because the trigger price and payoff strike differ, the terminal payoff is
# discontinuous at the trigger price.
#
# In this example the payoff jumps up from zero to approximately 20 when the
# terminal stock price crosses the trigger from below.
#
# The delta therefore develops an increasingly large positive spike around
# the trigger price as expiry approaches.
# ============================================================================

print("\n" + LINE)
print("5. GAP CALL DELTA NEAR EXPIRY")
print(LINE)

plt.figure()

for expiry_dt_near, maturity_label in expiry_dates:

    gap_call = EquityGapOption(
        expiry_dt_near,
        trigger_price,
        payoff_strike,
        OptionTypes.EUROPEAN_CALL,
    )

    delta_values = []

    for stock_price in stock_grid:

        delta = gap_call.delta(
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

    plt.plot(
        stock_grid,
        delta_values,
        label=maturity_label,
    )

plt.axvline(
    trigger_price,
    linestyle="--",
    linewidth=1.0,
    label=f"Trigger Price = {trigger_price:.0f}",
)

plt.axvline(
    payoff_strike,
    linestyle=":",
    linewidth=1.0,
    label=f"Payoff Strike = {payoff_strike:.0f}",
)

plt.axhline(
    0.0,
    linestyle=":",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Delta")

plt.title(
    "Gap Call Delta Near Expiry"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. ANALYTIC VERSUS MONTE CARLO
# ============================================================================
# What this section demonstrates:
#
# Compares the analytic and Monte Carlo values of a European equity gap call
# over a range of stock prices.
#
# The trigger price determines whether the payoff is activated and the payoff
# strike determines the amount paid.
#
# The Monte Carlo values should fluctuate around the analytic values because
# of sampling error.
# ============================================================================

print("\n" + LINE)
print("6. ANALYTIC VERSUS MONTE CARLO")
print(LINE)

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

trigger_price = 100.0
payoff_strike = 80.0

gap_call = EquityGapOption(
    expiry_dt,
    trigger_price,
    payoff_strike,
    OptionTypes.EUROPEAN_CALL,
)

stock_grid_mc = np.linspace(
    70.0,
    130.0,
    13,
)

num_paths = 100000
seed = 42

analytic_values = []
mc_values = []

print(
    f"{'S':>10}"
    f"{'ANALYTIC':>16}"
    f"{'MONTE CARLO':>16}"
    f"{'DIFFERENCE':>16}"
)

print("-" * 58)

for stock_price in stock_grid_mc:

    analytic_value = gap_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    mc_value = gap_call.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
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
    stock_grid_mc,
    analytic_values,
    marker="o",
    label="Analytic",
)

plt.plot(
    stock_grid_mc,
    mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.axvline(
    trigger_price,
    linestyle="--",
    linewidth=1.0,
    label=f"Trigger Price = {trigger_price:.0f}",
)

plt.axvline(
    payoff_strike,
    linestyle=":",
    linewidth=1.0,
    label=f"Payoff Strike = {payoff_strike:.0f}",
)

plt.axhline(
    0.0,
    linestyle="--",
    linewidth=1.0,
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")

plt.title(
    "Gap Call Analytic and Monte Carlo Values"
)

plt.grid(True)
plt.legend()
plt.show()
