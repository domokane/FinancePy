# ============================================================================
# FINANCEPY EXAMPLES - EquityCliquetOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# Demonstrates:
#   * call cliquets only
#   * PRICE and RETURN payoff types
#   * PERIODIC and MATURITY payment timing
#   * independent Monte Carlo checks of the analytic valuation
#
# The Monte Carlo check assumes valuation on the cliquet start/reset date.
# This avoids requiring historical reset fixings or already-realised payoffs.
# ============================================================================

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes
from financepy.utils.global_types import PaymentTimingTypes
from financepy.utils.global_types import CliquetTypes
from financepy.utils.frequency import FrequencyTypes

from financepy.products.equity.equity_cliquet_option import EquityCliquetOption
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve


# ============================================================================
# COMMON MARKET DATA
# ============================================================================

# Use value_dt == start_dt so every payoff/payment convention can be checked
# without supplying historical reset fixings.
value_dt = Date(1, 1, 2015)
start_dt = value_dt
final_expiry_dt = Date(1, 1, 2017)

stock_price = 100.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02
freq_type = FrequencyTypes.QUARTERLY
notional = 1.0

model = BlackScholes(volatility)
discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)

option_types = [
    OptionTypes.EUROPEAN_CALL,
]

payoff_types = [
    CliquetTypes.PRICE,
    CliquetTypes.RETURN,
]

payment_timings = [
    PaymentTimingTypes.PERIODIC,
    PaymentTimingTypes.MATURITY,
]


def short_name(x):
    """Return the enum member name without its class prefix."""
    return getattr(x, "name", str(x).split(".")[-1])


# ============================================================================
# 1. ALL PAYOFF / PAYMENT-TIMING COMBINATIONS
# ============================================================================

print("\n" + "=" * 108)
print("1. ANALYTIC VALUES - ALL CLIQUET CONFIGURATIONS")
print("=" * 108)

print(
    f"{'OPTION':<12}"
    f"{'PAYOFF':<12}"
    f"{'PAYMENT':<12}"
    f"{'VALUE':>18}"
)
print("-" * 54)

for opt_type in option_types:
    for payoff_type in payoff_types:
        for payoff_timing in payment_timings:

            cliquet = EquityCliquetOption(
                start_dt,
                final_expiry_dt,
                opt_type,
                freq_type,
                payoff_type,
                payoff_timing,
                notional,
            )

            value = cliquet.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            print(
                f"{short_name(opt_type):<12}"
                f"{short_name(payoff_type):<12}"
                f"{short_name(payoff_timing):<12}"
                f"{value:18.8f}"
            )

# ============================================================================
# 2. MONTE CARLO VALIDATION OF ALL 4 CALL CONFIGURATIONS
# ============================================================================

print("\n" + "=" * 108)
print("2. CALL CLIQUET: ANALYTIC VERSUS MONTE CARLO")
print("=" * 108)

num_paths = 500_000
seed = 4242

print(
    f"{'OPTION':<12}"
    f"{'PAYOFF':<12}"
    f"{'PAYMENT':<12}"
    f"{'ANALYTIC':>14}"
    f"{'MC':>14}"
    f"{'STD ERR':>12}"
    f"{'DIFF':>14}"
    f"{'Z':>10}"
)
print("-" * 110)

mc_results = []

for opt_type in option_types:
    for payoff_type in payoff_types:
        for payoff_timing in payment_timings:

            cliquet = EquityCliquetOption(
                start_dt,
                final_expiry_dt,
                opt_type,
                freq_type,
                payoff_type,
                payoff_timing,
                notional,
            )

            analytic = cliquet.value(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            mc_result = cliquet.value_mc(
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
                num_paths=num_paths,
                seed=seed,
            )

            mc_value = mc_result.value
            mc_stderr = mc_result.std_err

            diff = mc_value - analytic
            z_score = diff / mc_stderr if mc_stderr > 0.0 else np.nan
            passed = abs(z_score) <= 3.0

            mc_results.append(
                (
                    opt_type,
                    payoff_type,
                    payoff_timing,
                    analytic,
                    mc_value,
                    mc_stderr,
                    diff,
                    z_score,
                    passed,
                )
            )

            print(
                f"{short_name(opt_type):<12}"
                f"{short_name(payoff_type):<12}"
                f"{short_name(payoff_timing):<12}"
                f"{analytic:14.8f}"
                f"{mc_value:14.8f}"
                f"{mc_stderr:12.8f}"
                f"{diff:14.8f}"
                f"{z_score:10.3f}"
            )

print("\n3-sigma MC checks:")
for result in mc_results:
    opt_type, payoff_type, payoff_timing, _, _, _, _, z_score, passed = result
    status = "PASS" if passed else "FAIL"
    print(
        f"  {short_name(opt_type):<12} "
        f"{short_name(payoff_type):<8} "
        f"{short_name(payoff_timing):<10} "
        f"Z={z_score:8.3f}  {status}"
    )


# ============================================================================
# 3. VALUE VERSUS STOCK PRICE
# ============================================================================

print("\n" + "=" * 108)
print("3. CLIQUET VALUE VERSUS STOCK PRICE")
print("=" * 108)

stock_prices = np.linspace(50.0, 150.0, 21)

for payoff_type in payoff_types:

    plt.figure()

    for payoff_timing in payment_timings:

        cliquet = EquityCliquetOption(
            start_dt,
            final_expiry_dt,
            OptionTypes.EUROPEAN_CALL,
            freq_type,
            payoff_type,
            payoff_timing,
            notional,
        )

        values = []

        for s in stock_prices:
            values.append(
                cliquet.value(
                    value_dt,
                    s,
                    discount_curve,
                    dividend_curve,
                    model,
                )
            )

        plt.plot(
            stock_prices,
            values,
            marker="o",
            label=short_name(payoff_timing),
        )

    plt.xlabel("Stock Price")
    plt.ylabel("Cliquet Option Value")
    plt.title(
        "Call Cliquet Value versus Stock Price - "
        f"{short_name(payoff_type)}"
    )
    plt.grid(True)
    plt.legend(title="Payment Timing")
    plt.show()


# ============================================================================
# 4. VALUE VERSUS VOLATILITY
# ============================================================================

print("\n" + "=" * 108)
print("4. CLIQUET VALUE VERSUS VOLATILITY")
print("=" * 108)

volatilities = np.linspace(0.05, 0.50, 10)

for payoff_type in payoff_types:

    plt.figure()

    for payoff_timing in payment_timings:

        cliquet = EquityCliquetOption(
            start_dt,
            final_expiry_dt,
            OptionTypes.EUROPEAN_CALL,
            freq_type,
            payoff_type,
            payoff_timing,
            notional,
        )

        values = []

        for vol in volatilities:
            test_model = BlackScholes(vol)
            values.append(
                cliquet.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    test_model,
                )
            )

        plt.plot(
            volatilities,
            values,
            marker="o",
            label=short_name(payoff_timing),
        )

    plt.xlabel("Volatility")
    plt.ylabel("Cliquet Option Value")
    plt.title(
        "Call Cliquet Value versus Volatility - "
        f"{short_name(payoff_type)}"
    )
    plt.grid(True)
    plt.legend(title="Payment Timing")
    plt.show()


# ============================================================================
# 5. VALUE VERSUS INTEREST RATE
# ============================================================================

print("\n" + "=" * 108)
print("5. CLIQUET VALUE VERSUS INTEREST RATE")
print("=" * 108)

interest_rates = np.linspace(0.00, 0.10, 11)

for payoff_type in payoff_types:

    plt.figure()

    for payoff_timing in payment_timings:

        cliquet = EquityCliquetOption(
            start_dt,
            final_expiry_dt,
            OptionTypes.EUROPEAN_CALL,
            freq_type,
            payoff_type,
            payoff_timing,
            notional,
        )

        values = []

        for rate in interest_rates:
            test_discount_curve = FlatDiscountCurve(value_dt, rate)
            values.append(
                cliquet.value(
                    value_dt,
                    stock_price,
                    test_discount_curve,
                    dividend_curve,
                    model,
                )
            )

        plt.plot(
            interest_rates,
            values,
            marker="o",
            label=short_name(payoff_timing),
        )

    plt.xlabel("Interest Rate")
    plt.ylabel("Cliquet Option Value")
    plt.title(
        "Call Cliquet Value versus Interest Rate - "
        f"{short_name(payoff_type)}"
    )
    plt.grid(True)
    plt.legend(title="Payment Timing")
    plt.show()


# ============================================================================
# 6. VALUE VERSUS DIVIDEND YIELD
# ============================================================================

print("\n" + "=" * 108)
print("6. CLIQUET VALUE VERSUS DIVIDEND YIELD")
print("=" * 108)

dividend_yields = np.linspace(0.00, 0.10, 11)

for payoff_type in payoff_types:

    plt.figure()

    for payoff_timing in payment_timings:

        cliquet = EquityCliquetOption(
            start_dt,
            final_expiry_dt,
            OptionTypes.EUROPEAN_CALL,
            freq_type,
            payoff_type,
            payoff_timing,
            notional,
        )

        values = []

        for div_yield in dividend_yields:
            test_dividend_curve = FlatDiscountCurve(value_dt, div_yield)
            values.append(
                cliquet.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    test_dividend_curve,
                    model,
                )
            )

        plt.plot(
            dividend_yields,
            values,
            marker="o",
            label=short_name(payoff_timing),
        )

    plt.xlabel("Dividend Yield")
    plt.ylabel("Cliquet Option Value")
    plt.title(
        "Call Cliquet Value versus Dividend Yield - "
        f"{short_name(payoff_type)}"
    )
    plt.grid(True)
    plt.legend(title="Payment Timing")
    plt.show()


# ============================================================================
# 7. VALUE VERSUS RESET FREQUENCY
# ============================================================================

print("\n" + "=" * 108)
print("7. CLIQUET VALUE VERSUS RESET FREQUENCY")
print("=" * 108)

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

x = np.arange(len(frequencies))

for payoff_type in payoff_types:

    plt.figure()

    for payoff_timing in payment_timings:

        values = []

        for frequency in frequencies:
            cliquet = EquityCliquetOption(
                start_dt,
                final_expiry_dt,
                OptionTypes.EUROPEAN_CALL,
                frequency,
                payoff_type,
                payoff_timing,
                notional,
            )

            values.append(
                cliquet.value(
                    value_dt,
                    stock_price,
                    discount_curve,
                    dividend_curve,
                    model,
                )
            )

        plt.plot(
            x,
            values,
            marker="o",
            label=short_name(payoff_timing),
        )

    plt.xticks(x, frequency_labels)
    plt.xlabel("Reset Frequency")
    plt.ylabel("Cliquet Option Value")
    plt.title(
        "Call Cliquet Value versus Reset Frequency - "
        f"{short_name(payoff_type)}"
    )
    plt.grid(True)
    plt.legend(title="Payment Timing")
    plt.show()
