# ============================================================================
# FINANCEPY EXAMPLES - EquitySwap
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

"""
FinancePy example: EquitySwap

This is an executable example

It demonstrates how to build EquitySwap and EquitySwapLeg objects and then
uses them in a collection of exploratory examples:

    1. Create and print an equity swap.
    2. Value a clean swap at inception.
    3. Revalue a swap after inception using a previous fixing.
    4. Value an equity leg with and without dividends.
    5. Sweep the current stock price and plot the swap value.
    6. Sweep dividend yield and plot the swap value.
    7. Sweep the floating/index rate and plot the swap value.
    8. Sweep the floating-leg spread and plot the swap value.
    9. Compare receive-equity and pay-equity structures.
   10. Display receive/pay value relationships.
   11. Plot stock-price sensitivity for several dividend yields.
   12. Run lightweight numerical sanity checks.

The checks at the end are deliberately part of the example: they make it easy
to use this file interactively while developing or changing EquitySwap code.
They are not intended to replace FinancePy's formal unit tests.
"""

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.global_types import SwapTypes, ReturnTypes
from financepy.utils.global_vars import ONE_MILLION

from financepy.products.equity.equity_swap import EquitySwap
from financepy.products.equity.equity_swap_leg import EquitySwapLeg

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


# ============================================================================
# HELPERS
# ============================================================================

def print_header(number, title):
    print()
    print("=" * 78)
    print(f"{number}. {title}")
    print("=" * 78)


def make_equity_swap(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=130.0,
    notional=ONE_MILLION,
    rate_spread=0.0,
    return_type=ReturnTypes.TOTAL_RETURN,
    freq_type=FrequencyTypes.SEMI_ANNUAL,
    dc_type=DayCountTypes.THIRTY_360_BOND,
    payment_lag=0,
):
    """Construct a standard equity swap used throughout this example."""

    stock_qty = notional / stock_strike

    return EquitySwap(
        effective_dt,
        maturity_dt,
        leg_type,
        freq_type,
        dc_type,
        stock_strike,
        stock_qty,
        payment_lag,
        return_type,
        freq_type,
        dc_type,
        rate_spread,
        payment_lag,
        CalendarTypes.TARGET,
        BusDayAdjustTypes.FOLLOWING,
        DateGenRuleTypes.BACKWARD,
    )


def make_equity_leg(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=130.0,
    notional=ONE_MILLION,
    return_type=ReturnTypes.TOTAL_RETURN,
    freq_type=FrequencyTypes.SEMI_ANNUAL,
    dc_type=DayCountTypes.THIRTY_360_BOND,
    payment_lag=0,
):
    """Construct an equity swap leg."""

    stock_qty = notional / stock_strike

    return EquitySwapLeg(
        effective_dt,
        maturity_dt,
        leg_type,
        freq_type,
        dc_type,
        stock_strike,
        stock_qty,
        payment_lag,
        return_type,
        CalendarTypes.TARGET,
        BusDayAdjustTypes.FOLLOWING,
        DateGenRuleTypes.BACKWARD,
    )


# ============================================================================
# COMMON MARKET / CONTRACT DATA
# ============================================================================

effective_dt = Date(13, 2, 2018)
maturity_dt = effective_dt.add_months(24)

stock_strike = 130.0
notional = ONE_MILLION

discount_rate = 0.05
index_rate = 0.03
dividend_rate = 0.02

discount_curve = FlatDiscountCurve(
    effective_dt,
    discount_rate,
)

index_curve = FlatDiscountCurve(
    effective_dt,
    index_rate,
)

dividend_curve = FlatDiscountCurve(
    effective_dt,
    dividend_rate,
)


# ============================================================================
# 1. CREATE AN EQUITY SWAP
# ============================================================================

print_header(1, "CREATE AN EQUITY SWAP")

equity_swap = make_equity_swap(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=stock_strike,
    notional=notional,
    return_type=ReturnTypes.TOTAL_RETURN,
)

print(equity_swap)


# ============================================================================
# 2. VALUE AT INCEPTION
# ============================================================================

print_header(2, "VALUE AT INCEPTION")

# For this clean inception test use the same curve for discounting and
# floating-rate projection and zero dividend yield, matching the original
# FinancePy regression example.

inception_discount_curve = FlatDiscountCurve(
    effective_dt,
    0.05,
)

inception_dividend_curve = FlatDiscountCurve(
    effective_dt,
    0.00,
)

inception_value = equity_swap.value(
    effective_dt,
    inception_discount_curve,
    inception_discount_curve,
    inception_dividend_curve,
)

print(f"Swap value : {inception_value:14.8f}")

print(
    "Near zero at inception :",
    abs(inception_value) < 1.0e-5,
)


# ============================================================================
# 3. REVALUE AFTER INCEPTION
# ============================================================================

print_header(3, "REVALUE AFTER INCEPTION")

effective_dt_1y = Date(13, 2, 2018)
value_dt_1y = effective_dt_1y.add_months(6)
maturity_dt_1y = effective_dt_1y.add_months(12)

stock_strike_1y = 125.0
notional_1y = ONE_MILLION

discount_rate_1y = 0.05

discount_curve_1y = FlatDiscountCurve(
    value_dt_1y,
    discount_rate_1y,
)

dividend_curve_1y = FlatDiscountCurve(
    value_dt_1y,
    0.0,
)

index_curve_1y = discount_curve_1y

# Floating rate fixed at the previous reset date.
index_curve_first = FlatDiscountCurve(
    effective_dt_1y,
    discount_rate_1y,
)

index_alpha_first = DayCount(
    index_curve_first.time_dc_type
).year_frac(
    effective_dt_1y,
    maturity_dt_1y,
)[0]

first_fixing = (
    index_curve_first.df(effective_dt_1y)
    / index_curve_first.df(maturity_dt_1y)
    - 1.0
) / index_alpha_first

# Forward rate from valuation date to maturity.
index_curve_period = FlatDiscountCurve(
    value_dt_1y,
    discount_rate_1y,
)

index_alpha_period = DayCount(
    index_curve_period.dc_type
).year_frac(
    value_dt_1y,
    maturity_dt_1y,
)[0]

period_fixing = (
    index_curve_period.df(value_dt_1y)
    / index_curve_period.df(maturity_dt_1y)
    - 1.0
) / index_alpha_period

# Stock price for which the equity and floating legs balance.
stock_price_zero = (
    stock_strike_1y
    * (1.0 + first_fixing * index_alpha_first)
    / (1.0 + period_fixing * index_alpha_period)
)

equity_swap_1y = make_equity_swap(
    effective_dt_1y,
    maturity_dt_1y,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=stock_strike_1y,
    notional=notional_1y,
    freq_type=FrequencyTypes.ANNUAL,
)

value_zero = equity_swap_1y.value(
    value_dt_1y,
    discount_curve_1y,
    index_curve_1y,
    dividend_curve_1y,
    stock_price_zero,
    first_fixing,
)

print(f"First fixing       : {first_fixing:14.8f}")
print(f"Period fixing      : {period_fixing:14.8f}")
print(f"Zero-value stock   : {stock_price_zero:14.8f}")
print(f"Swap value         : {value_zero:14.8f}")


# ============================================================================
# 4. EQUITY LEG WITH DIVIDENDS
# ============================================================================

print_header(4, "EQUITY LEG WITH DIVIDENDS")

equity_leg = make_equity_leg(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=stock_strike,
    notional=notional,
    return_type=ReturnTypes.TOTAL_RETURN,
)

zero_dividend_curve = FlatDiscountCurve(
    effective_dt,
    0.0,
)

value_no_divs = equity_leg.value(
    effective_dt,
    discount_curve,
    index_curve,
    zero_dividend_curve,
    stock_strike,
)

value_with_divs = equity_leg.value(
    effective_dt,
    discount_curve,
    index_curve,
    dividend_curve,
    stock_strike,
)

print(f"Equity leg, q=0%   : {value_no_divs:14.8f}")
print(f"Equity leg, q=2%   : {value_with_divs:14.8f}")
print(f"Difference         : {value_with_divs - value_no_divs:14.8f}")


# ============================================================================
# 5. SWAP VALUE VERSUS STOCK PRICE
# ============================================================================

print_header(5, "SWAP VALUE VERSUS STOCK PRICE")

stock_prices = np.linspace(
    80.0,
    180.0,
    41,
)

stock_values = []

for stock_price in stock_prices:

    value = equity_swap.value(
        effective_dt,
        discount_curve,
        index_curve,
        dividend_curve,
        stock_price,
    )

    stock_values.append(value)

print(
    f"{'STOCK':>12}"
    f"{'SWAP VALUE':>18}"
)

print("-" * 30)

for stock_price, value in zip(
    stock_prices[::5],
    stock_values[::5],
):
    print(
        f"{stock_price:12.4f}"
        f"{value:18.4f}"
    )

fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    stock_values,
    marker="o",
    markersize=3,
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.axvline(
    stock_strike,
    linestyle=":",
    label="Initial stock price",
)

ax.set_xlabel("Current stock price")
ax.set_ylabel("Equity swap value")
ax.set_title("Equity Swap Value versus Stock Price")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 6. SWAP VALUE VERSUS DIVIDEND YIELD
# ============================================================================

print_header(6, "SWAP VALUE VERSUS DIVIDEND YIELD")

dividend_rates = np.linspace(
    0.0,
    0.08,
    17,
)

dividend_values = []

for q in dividend_rates:

    q_curve = FlatDiscountCurve(
        effective_dt,
        q,
    )

    value = equity_swap.value(
        effective_dt,
        discount_curve,
        index_curve,
        q_curve,
        stock_strike,
    )

    dividend_values.append(value)

print(
    f"{'DIV YIELD':>12}"
    f"{'SWAP VALUE':>18}"
)

print("-" * 30)

for q, value in zip(
    dividend_rates,
    dividend_values,
):
    print(
        f"{q:12.4%}"
        f"{value:18.4f}"
    )

fig, ax = plt.subplots()

ax.plot(
    100.0 * dividend_rates,
    dividend_values,
    marker="o",
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel("Dividend yield (%)")
ax.set_ylabel("Equity swap value")
ax.set_title("Equity Swap Value versus Dividend Yield")
ax.grid(True)

plt.show()


# ============================================================================
# 7. SWAP VALUE VERSUS INDEX RATE
# ============================================================================

print_header(7, "SWAP VALUE VERSUS INDEX RATE")

index_rates = np.linspace(
    0.0,
    0.08,
    17,
)

index_values = []

for rate in index_rates:

    test_index_curve = FlatDiscountCurve(
        effective_dt,
        rate,
    )

    value = equity_swap.value(
        effective_dt,
        discount_curve,
        test_index_curve,
        dividend_curve,
        stock_strike,
    )

    index_values.append(value)

print(
    f"{'INDEX RATE':>12}"
    f"{'SWAP VALUE':>18}"
)

print("-" * 30)

for rate, value in zip(
    index_rates,
    index_values,
):
    print(
        f"{rate:12.4%}"
        f"{value:18.4f}"
    )

fig, ax = plt.subplots()

ax.plot(
    100.0 * index_rates,
    index_values,
    marker="o",
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel("Index rate (%)")
ax.set_ylabel("Equity swap value")
ax.set_title("Equity Swap Value versus Index Rate")
ax.grid(True)

plt.show()


# ============================================================================
# 8. SWAP VALUE VERSUS RATE SPREAD
# ============================================================================

print_header(8, "SWAP VALUE VERSUS RATE SPREAD")

spreads = np.linspace(
    -0.02,
    0.02,
    17,
)

spread_values = []

for spread in spreads:

    swap_with_spread = make_equity_swap(
        effective_dt,
        maturity_dt,
        leg_type=SwapTypes.RECEIVE,
        stock_strike=stock_strike,
        notional=notional,
        rate_spread=spread,
        return_type=ReturnTypes.TOTAL_RETURN,
    )

    value = swap_with_spread.value(
        effective_dt,
        discount_curve,
        index_curve,
        dividend_curve,
        stock_strike,
    )

    spread_values.append(value)

print(
    f"{'SPREAD':>12}"
    f"{'SWAP VALUE':>18}"
)

print("-" * 30)

for spread, value in zip(
    spreads,
    spread_values,
):
    print(
        f"{spread:12.4%}"
        f"{value:18.4f}"
    )

fig, ax = plt.subplots()

ax.plot(
    10000.0 * spreads,
    spread_values,
    marker="o",
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel("Rate spread (bp)")
ax.set_ylabel("Equity swap value")
ax.set_title("Equity Swap Value versus Rate Spread")
ax.grid(True)

plt.show()


# ============================================================================
# 9. RECEIVE VERSUS PAY EQUITY
# ============================================================================

print_header(9, "RECEIVE VERSUS PAY EQUITY")

receive_swap = make_equity_swap(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.RECEIVE,
    stock_strike=stock_strike,
    notional=notional,
)

pay_swap = make_equity_swap(
    effective_dt,
    maturity_dt,
    leg_type=SwapTypes.PAY,
    stock_strike=stock_strike,
    notional=notional,
)

receive_values = []
pay_values = []

for stock_price in stock_prices:

    receive_value = receive_swap.value(
        effective_dt,
        discount_curve,
        index_curve,
        dividend_curve,
        stock_price,
    )

    pay_value = pay_swap.value(
        effective_dt,
        discount_curve,
        index_curve,
        dividend_curve,
        stock_price,
    )

    receive_values.append(receive_value)
    pay_values.append(pay_value)

fig, ax = plt.subplots()

ax.plot(
    stock_prices,
    receive_values,
    label="Receive equity",
)

ax.plot(
    stock_prices,
    pay_values,
    label="Pay equity",
)

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel("Current stock price")
ax.set_ylabel("Equity swap value")
ax.set_title("Receive versus Pay Equity Swap")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 10. RECEIVE/PAY SIGN CHECK
# ============================================================================

print_header(10, "RECEIVE/PAY SIGN CHECK")

print(
    f"{'STOCK':>12}"
    f"{'RECEIVE':>18}"
    f"{'PAY':>18}"
    f"{'SUM':>18}"
)

print("-" * 66)

for stock_price, receive_value, pay_value in zip(
    stock_prices[::5],
    receive_values[::5],
    pay_values[::5],
):

    value_sum = receive_value + pay_value

    print(
        f"{stock_price:12.4f}"
        f"{receive_value:18.4f}"
        f"{pay_value:18.4f}"
        f"{value_sum:18.8f}"
    )


# ============================================================================
# 11. TWO-DIMENSIONAL STOCK / DIVIDEND SENSITIVITY
# ============================================================================

print_header(
    11,
    "STOCK / DIVIDEND SENSITIVITY",
)

stock_grid = np.linspace(
    90.0,
    170.0,
    17,
)

dividend_grid = np.array(
    [0.00, 0.01, 0.02, 0.04, 0.06]
)

fig, ax = plt.subplots()

for q in dividend_grid:

    q_curve = FlatDiscountCurve(
        effective_dt,
        q,
    )

    values = []

    for stock_price in stock_grid:

        value = equity_swap.value(
            effective_dt,
            discount_curve,
            index_curve,
            q_curve,
            stock_price,
        )

        values.append(value)

    ax.plot(
        stock_grid,
        values,
        marker="o",
        label=f"q = {q:.0%}",
    )

ax.axhline(
    0.0,
    linestyle="--",
)

ax.set_xlabel("Current stock price")
ax.set_ylabel("Equity swap value")
ax.set_title("Stock Price and Dividend-Yield Sensitivity")
ax.grid(True)
ax.legend()

plt.show()


# ============================================================================
# 12. BASIC SANITY CHECKS
# ============================================================================

print_header(12, "BASIC SANITY CHECKS")

checks = {
    "Stock-price sweep finite": np.all(np.isfinite(stock_values)),
    "Dividend sweep finite": np.all(np.isfinite(dividend_values)),
    "Index-rate sweep finite": np.all(np.isfinite(index_values)),
    "Spread sweep finite": np.all(np.isfinite(spread_values)),
    "Receive-equity sweep finite": np.all(np.isfinite(receive_values)),
    "Pay-equity sweep finite": np.all(np.isfinite(pay_values)),
}

for name, passed in checks.items():
    print(f"{name:<36}: {passed}")

print()
print("Example completed.")
