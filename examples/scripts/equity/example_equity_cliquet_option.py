# ============================================================================
# FINANCEPY EXAMPLES - EquityCliquetOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

# Allow this example to run directly from its category folder.

import numpy as np
import matplotlib.pyplot as plt

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.utils.format_graphs import set_plot_style


set_plot_style()


# ============================================================================
# 1. EQUITY CLIQUET OPTION
# ============================================================================
# What this section demonstrates:
# Values a cliquet option using the supplied market data and model inputs.
#
# A cliquet consists of a sequence of forward-start options. At each reset
# date a new option effectively begins, with its strike determined by the
# stock level at that reset date.

print("\n" + "=" * 78)
print("1. EQUITY CLIQUET OPTION")
print("=" * 78)

start_dt = Date(1, 1, 2014)
final_expiry_dt = Date(1, 1, 2017)

freq_type = FrequencyTypes.QUARTERLY
opt_type = OptionTypes.EUROPEAN_CALL

cliquet_option = EquityCliquetOption(
    start_dt,
    final_expiry_dt,
    opt_type,
    freq_type,
)

value_dt = Date(1, 1, 2015)

stock_price = 100.0
volatility = 0.20
interest_rate = 0.05
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

value = cliquet_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
)

print("LABEL", "VALUE")
print("FINANCEPY", value)


# ============================================================================
# 2. CALL AND PUT CLIQUET OPTIONS
# ============================================================================
# What this section demonstrates:
# Compares otherwise identical call and put cliquet options.
#
# The comparison illustrates how the direction of the individual forward-
# starting option payoffs affects the total cliquet value.

print("\n" + "=" * 78)
print("2. CALL AND PUT CLIQUET OPTIONS")
print("=" * 78)

print(
    f"{'OPTION TYPE':<20}"
    f"{'VALUE':>16}"
)

print("-" * 36)

for opt_type in [
    OptionTypes.EUROPEAN_CALL,
    OptionTypes.EUROPEAN_PUT,
]:

    cliquet_option = EquityCliquetOption(
        start_dt,
        final_expiry_dt,
        opt_type,
        freq_type,
    )

    value = cliquet_option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )

    print(
        f"{str(opt_type):<20}"
        f"{value:16.8f}"
    )


# ============================================================================
# 3. CLIQUET VALUE VERSUS STOCK PRICE
# ============================================================================
# What this section demonstrates:
# Shows how the cliquet value changes with the current stock price.
#
# This is particularly useful for a cliquet because its future strikes are
# reset relative to future stock levels rather than being fixed today.

print("\n" + "=" * 78)
print("3. CLIQUET VALUE VERSUS STOCK PRICE")
print("=" * 78)

stock_prices = np.linspace(
    50.0,
    150.0,
    21,
)

call_option = EquityCliquetOption(
    start_dt,
    final_expiry_dt,
    OptionTypes.EUROPEAN_CALL,
    freq_type,
)

put_option = EquityCliquetOption(
    start_dt,
    final_expiry_dt,
    OptionTypes.EUROPEAN_PUT,
    freq_type,
)

call_values = []
put_values = []

print(
    f"{'STOCK':>10}"
    f"{'CALL':>16}"
    f"{'PUT':>16}"
)

print("-" * 42)

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

    print(
        f"{stock_price:10.2f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

call_values = np.asarray(call_values)
put_values = np.asarray(put_values)

plt.figure()

plt.plot(
    stock_prices,
    call_values,
    marker="o",
    label="Call Cliquet",
)

plt.plot(
    stock_prices,
    put_values,
    marker="o",
    label="Put Cliquet",
)

plt.xlabel("Stock Price")
plt.ylabel("Cliquet Option Value")
plt.title("Cliquet Option Value versus Stock Price")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. CLIQUET VALUE VERSUS VOLATILITY
# ============================================================================
# What this section demonstrates:
# Shows the sensitivity of the cliquet value to volatility.
#
# Since the cliquet contains a sequence of option-like payoffs, volatility is
# an important determinant of its value.

print("\n" + "=" * 78)
print("4. CLIQUET VALUE VERSUS VOLATILITY")
print("=" * 78)

stock_price = 100.0

volatilities = np.linspace(
    0.05,
    0.50,
    10,
)

call_values = []
put_values = []

print(
    f"{'VOLATILITY':>12}"
    f"{'CALL':>16}"
    f"{'PUT':>16}"
)

print("-" * 44)

for volatility in volatilities:

    model = BlackScholes(
        volatility,
    )

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

    print(
        f"{volatility:12.4f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

plt.figure()

plt.plot(
    volatilities,
    call_values,
    marker="o",
    label="Call Cliquet",
)

plt.plot(
    volatilities,
    put_values,
    marker="o",
    label="Put Cliquet",
)

plt.xlabel("Volatility")
plt.ylabel("Cliquet Option Value")
plt.title("Cliquet Option Value versus Volatility")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. CLIQUET VALUE VERSUS INTEREST RATE
# ============================================================================
# What this section demonstrates:
# Shows how the cliquet value responds to changes in the interest-rate curve
# while holding the other market inputs fixed.

print("\n" + "=" * 78)
print("5. CLIQUET VALUE VERSUS INTEREST RATE")
print("=" * 78)

volatility = 0.20
model = BlackScholes(volatility)

interest_rates = np.linspace(
    0.00,
    0.10,
    11,
)

call_values = []
put_values = []

print(
    f"{'RATE':>12}"
    f"{'CALL':>16}"
    f"{'PUT':>16}"
)

print("-" * 44)

for interest_rate in interest_rates:

    test_discount_curve = FlatDiscountCurve(
        value_dt,
        interest_rate,
    )

    call_value = call_option.value(
        value_dt,
        stock_price,
        test_discount_curve,
        dividend_curve,
        model,
    )

    put_value = put_option.value(
        value_dt,
        stock_price,
        test_discount_curve,
        dividend_curve,
        model,
    )

    call_values.append(call_value)
    put_values.append(put_value)

    print(
        f"{interest_rate:12.4f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

plt.figure()

plt.plot(
    interest_rates,
    call_values,
    marker="o",
    label="Call Cliquet",
)

plt.plot(
    interest_rates,
    put_values,
    marker="o",
    label="Put Cliquet",
)

plt.xlabel("Interest Rate")
plt.ylabel("Cliquet Option Value")
plt.title("Cliquet Option Value versus Interest Rate")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. CLIQUET VALUE VERSUS DIVIDEND YIELD
# ============================================================================
# What this section demonstrates:
# Shows the effect of changing the dividend yield while keeping the other
# market inputs fixed.
#
# Dividend yield changes the risk-neutral growth rate of the stock and
# therefore affects the forward-start option components of the cliquet.

print("\n" + "=" * 78)
print("6. CLIQUET VALUE VERSUS DIVIDEND YIELD")
print("=" * 78)

interest_rate = 0.05

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_yields = np.linspace(
    0.00,
    0.10,
    11,
)

call_values = []
put_values = []

print(
    f"{'DIV YIELD':>12}"
    f"{'CALL':>16}"
    f"{'PUT':>16}"
)

print("-" * 44)

for dividend_yield in dividend_yields:

    test_dividend_curve = FlatDiscountCurve(
        value_dt,
        dividend_yield,
    )

    call_value = call_option.value(
        value_dt,
        stock_price,
        discount_curve,
        test_dividend_curve,
        model,
    )

    put_value = put_option.value(
        value_dt,
        stock_price,
        discount_curve,
        test_dividend_curve,
        model,
    )

    call_values.append(call_value)
    put_values.append(put_value)

    print(
        f"{dividend_yield:12.4f}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

plt.figure()

plt.plot(
    dividend_yields,
    call_values,
    marker="o",
    label="Call Cliquet",
)

plt.plot(
    dividend_yields,
    put_values,
    marker="o",
    label="Put Cliquet",
)

plt.xlabel("Dividend Yield")
plt.ylabel("Cliquet Option Value")
plt.title("Cliquet Option Value versus Dividend Yield")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. CLIQUET VALUE VERSUS RESET FREQUENCY
# ============================================================================
# What this section demonstrates:
# Compares cliquet values for different reset frequencies.
#
# Changing the frequency changes both the number and length of the
# forward-start option periods making up the cliquet.

print("\n" + "=" * 78)
print("7. CLIQUET VALUE VERSUS RESET FREQUENCY")
print("=" * 78)

dividend_yield = 0.02

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

frequencies = [
    FrequencyTypes.ANNUAL,
    FrequencyTypes.SEMI_ANNUAL,
    FrequencyTypes.QUARTERLY,
    FrequencyTypes.MONTHLY,
]

frequency_labels = [
    "Annual",
    "Semi-Annual",
    "Quarterly",
    "Monthly",
]

call_values = []
put_values = []

print(
    f"{'FREQUENCY':<16}"
    f"{'CALL':>16}"
    f"{'PUT':>16}"
)

print("-" * 48)

for frequency, label in zip(
    frequencies,
    frequency_labels,
):

    call_option = EquityCliquetOption(
        start_dt,
        final_expiry_dt,
        OptionTypes.EUROPEAN_CALL,
        frequency,
    )

    put_option = EquityCliquetOption(
        start_dt,
        final_expiry_dt,
        OptionTypes.EUROPEAN_PUT,
        frequency,
    )

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

    print(
        f"{label:<16}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )

x = np.arange(
    len(frequencies),
)

plt.figure()

plt.plot(
    x,
    call_values,
    marker="o",
    label="Call Cliquet",
)

plt.plot(
    x,
    put_values,
    marker="o",
    label="Put Cliquet",
)

plt.xticks(
    x,
    frequency_labels,
)

plt.xlabel("Reset Frequency")
plt.ylabel("Cliquet Option Value")
plt.title("Cliquet Option Value versus Reset Frequency")
plt.grid(True)
plt.legend()
plt.show()
