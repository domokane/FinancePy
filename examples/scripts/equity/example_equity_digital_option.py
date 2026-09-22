# ============================================================================
# FINANCEPY EXAMPLES - EquityDigitalOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_digital_option import EquityDigitalOption
from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style
from financepy.utils.global_types import DigitalOptionTypes
from financepy.utils.global_types import OptionTypes

set_plot_style()


# ============================================================================
# 1. CASH-OR-NOTHING DIGITAL OPTION
# ============================================================================
# What this section demonstrates:
# Compares the Black-Scholes analytic value with Monte Carlo valuation for a
# cash-or-nothing digital call.
# Shows Monte Carlo convergence as the number of paths increases.

print("\n" + "=" * 78)
print("1. CASH-OR-NOTHING DIGITAL OPTION")
print("=" * 78)

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

stock_price = 100.0
strike_price = 100.0
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01

underlying_type = DigitalOptionTypes.CASH_OR_NOTHING

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)
model = BlackScholes(volatility)

call_option = EquityDigitalOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    underlying_type,
)

value_bs = call_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

num_paths_list = [
    1_000,
    2_000,
    5_000,
    10_000,
    20_000,
    50_000,
    100_000,
]

mc_values = []

print()
print("NUM PATHS       ANALYTIC             MC          ERROR       TIME")
print("-" * 70)

for num_paths in num_paths_list:

    start = time.time()

    value_mc = call_option.value_mc(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
    )

    duration = time.time() - start
    error = value_mc - value_bs

    mc_values.append(value_mc)

    print(
        f"{num_paths:9d}"
        f"{value_bs:15.8f}"
        f"{value_mc:15.8f}"
        f"{error:15.8f}"
        f"{duration:12.4f}"
    )


# ============================================================================
# 2. MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
# Shows how the Monte Carlo estimate fluctuates around the analytic
# Black-Scholes value as the number of simulated paths increases.

print("\n" + "=" * 78)
print("2. MONTE CARLO CONVERGENCE")
print("=" * 78)

fig, ax = plt.subplots()

ax.plot(
    num_paths_list,
    mc_values,
    marker="o",
    label="Monte Carlo",
)

ax.axhline(
    value_bs,
    linestyle="--",
    label="Analytic",
)

ax.set_xscale("log")
ax.set_xlabel("Number of Monte Carlo Paths")
ax.set_ylabel("Digital Option Value")
ax.set_title("Digital Option Monte Carlo Convergence")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 3. DIGITAL OPTION VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows the characteristic value profile of cash-or-nothing digital call
# and put options as the underlying stock price moves through the strike.

print("\n" + "=" * 78)
print("3. DIGITAL OPTION VALUE VERSUS STOCK PRICE")
print("=" * 78)

stock_prices = np.linspace(50.0, 150.0, 201)

call_option = EquityDigitalOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    underlying_type,
)

put_option = EquityDigitalOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_PUT,
    underlying_type,
)

call_values = []
put_values = []

for stock_price in stock_prices:

    call_value = call_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    put_value = put_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    call_values.append(call_value)
    put_values.append(put_value)

fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    call_values,
    label="Digital Call",
)

ax.plot(
    stock_prices,
    put_values,
    label="Digital Put",
)

ax.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Digital Option Value")
ax.set_title("Cash-or-Nothing Digital Option Value")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 4. DIGITAL OPTION GREEKS
# ============================================================================
# What this section demonstrates:
# Shows how delta, vega and theta behave around the strike. Digital options
# have particularly concentrated sensitivities near the strike because of
# their discontinuous terminal payoff.

print("\n" + "=" * 78)
print("4. DIGITAL OPTION GREEKS")
print("=" * 78)

stock_prices = np.linspace(60.0, 140.0, 161)

call_deltas = []
call_vegas = []
call_thetas = []

put_deltas = []
put_vegas = []
put_thetas = []

for stock_price in stock_prices:

    call_deltas.append(
        call_option.delta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

    call_vegas.append(
        call_option.vega(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

    call_thetas.append(
        call_option.theta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

    put_deltas.append(
        put_option.delta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

    put_vegas.append(
        put_option.vega(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )

    put_thetas.append(
        put_option.theta(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            model,
        )
    )


fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    call_deltas,
    label="Call Delta",
)

ax.plot(
    stock_prices,
    put_deltas,
    label="Put Delta",
)

ax.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Delta")
ax.set_title("Cash-or-Nothing Digital Option Delta")
ax.grid(True)
ax.legend()

plt.show()


fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    call_vegas,
    label="Call Vega",
)

ax.plot(
    stock_prices,
    put_vegas,
    label="Put Vega",
)

ax.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Vega")
ax.set_title("Cash-or-Nothing Digital Option Vega")
ax.grid(True)
ax.legend()

plt.show()


fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    call_thetas,
    label="Call Theta",
)

ax.plot(
    stock_prices,
    put_thetas,
    label="Put Theta",
)

ax.axvline(
    strike_price,
    linestyle="--",
    label="Strike",
)

ax.set_xlabel("Stock Price")
ax.set_ylabel("Theta")
ax.set_title("Cash-or-Nothing Digital Option Theta")
ax.grid(True)
ax.legend()

plt.show()
