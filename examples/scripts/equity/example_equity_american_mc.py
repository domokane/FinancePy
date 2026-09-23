# ============================================================================
# FINANCEPY EXAMPLES - EquityAmericanOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates:
#
#   1. European and American call/put valuation
#   2. Early-exercise premium
#   3. Comparison of FinancePy American-option pricing methods
#   4. CRR tree versus LSMC
#   5. Analytic/semi-analytic approximation comparison
#   6. Option value versus stock price
#   7. Early-exercise premium versus stock price
#   8. Approximation error relative to a CRR benchmark
#   9. Effect of dividends on American calls
#  10. Option value versus volatility
#  11. Option value as expiry approaches
#  12. Delta versus stock price
#  13. Gamma versus stock price
#  14. CRR tree convergence
#  15. CRR versus LSMC
#  16. Computational speed
#  17. Delta bump-and-revalue verification
#  18. Gamma bump-and-revalue verification
#
# The high-step CRR tree is used as a numerical benchmark for American
# exercise. Analytic/semi-analytic methods are compared with this benchmark
# where they are supported by the installed FinancePy version.
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes, BlackScholesTypes
from financepy.models.black_scholes_analytic import value as bs_value
from financepy.models.equity_crr_tree import crr_tree_val
from financepy.models.equity_lsmc import equity_lsmc, BoundaryFitTypes
from financepy.products.equity.equity_american_option import (
    EquityAmericanOption,
)


LINE = "=" * 100
SUBLINE = "-" * 100


# ============================================================================
# SUPPORTING FUNCTIONS
# ============================================================================


def option_value(
    option,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
):
    """Return the scalar option value."""

    return option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
    )


# ============================================================================


def build_curves(
    value_dt,
    interest_rate,
    dividend_yield,
):
    """Build flat interest-rate and dividend-yield curves."""

    discount_curve = FlatDiscountCurve(
        value_dt,
        interest_rate,
    )

    dividend_curve = FlatDiscountCurve(
        value_dt,
        dividend_yield,
    )

    return discount_curve, dividend_curve


# ============================================================================


def model_type_name(model_type):
    """Return a readable name for a FinancePy model type."""

    try:
        return model_type.name
    except AttributeError:
        return str(model_type)


# ============================================================================


def available_american_models(
    option,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    volatility,
):
    """
    Determine which BlackScholesTypes work through EquityAmericanOption.

    CRR_TREE and LSMC are handled separately because they require numerical
    configuration such as tree steps and Monte Carlo paths.
    """

    available = []

    for model_type in BlackScholesTypes:

        if model_type in (
            BlackScholesTypes.CRR_TREE,
            BlackScholesTypes.LSMC,
        ):
            continue

        try:

            model = BlackScholes(
                volatility,
                model_type,
            )

            value = option_value(
                option,
                value_dt,
                stock_price,
                discount_curve,
                dividend_curve,
                model,
            )

            if np.isfinite(value):

                available.append(
                    (
                        model_type,
                        value,
                    )
                )

        except Exception:
            pass

    return available


# ============================================================================


def replicate_ls_paper():
    """
    Reproduce the Longstaff-Schwartz comparison used in the original example.

    The LSMC value is compared with a CRR tree and the corresponding European
    Black-Scholes value.
    """

    amer_option_put_value = OptionTypes.AMERICAN_PUT.value

    opt_type_values = [
        amer_option_put_value,
    ]

    stock_prices = [
        36.0,
        38.0,
        40.0,
        42.0,
        44.0,
    ]

    volatilities = [
        0.20,
        0.40,
    ]

    times_to_expiry = [
        1.0,
        2.0,
    ]

    num_paths = 50000
    num_steps_per_year = 50

    r = 0.06
    q = 0.0
    k = 40.0

    poly_deg = 3
    fit_type_value = BoundaryFitTypes.HERMITE_E.value
    use_sobol = False
    seed = 1912

    print(
        f"{'S':>8}"
        f"{'VOL':>10}"
        f"{'T':>10}"
        f"{'TREE':>14}"
        f"{'EUROPEAN':>14}"
        f"{'LSMC':>14}"
    )

    print(SUBLINE)

    for opt_type_value in opt_type_values:

        for s in stock_prices:

            for v in volatilities:

                for t in times_to_expiry:

                    v_ls = equity_lsmc(
                        s,
                        r,
                        q,
                        v,
                        num_paths,
                        num_steps_per_year,
                        t,
                        opt_type_value,
                        k,
                        poly_deg,
                        fit_type_value,
                        use_sobol,
                        seed,
                    )

                    v_tree = crr_tree_val(
                        s,
                        r,
                        q,
                        v,
                        2000,
                        t,
                        opt_type_value,
                        k,
                        True,
                    )[0]

                    v_eur = bs_value(
                        s,
                        t,
                        k,
                        r,
                        q,
                        v,
                        OptionTypes.EUROPEAN_PUT.value,
                    )

                    print(
                        f"{s:8.2f}"
                        f"{v:10.4f}"
                        f"{t:10.4f}"
                        f"{v_tree:14.6f}"
                        f"{v_eur:14.6f}"
                        f"{v_ls:14.6f}"
                    )


# ============================================================================
# MARKET INPUTS
# ============================================================================

value_dt = Date(
    1,
    1,
    2016,
)

expiry_dt = Date(
    1,
    1,
    2017,
)

stock_price = 50.0
strike_price = 50.0

interest_rate = 0.06
dividend_yield = 0.04
volatility = 0.40

discount_curve, dividend_curve = build_curves(
    value_dt,
    interest_rate,
    dividend_yield,
)

num_steps = 500

crr_model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    num_steps,
)

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


# ============================================================================
# 1. EUROPEAN AND AMERICAN OPTION VALUATION
# ============================================================================

print("\n" + LINE)
print("1. EUROPEAN AND AMERICAN OPTION VALUATION")
print(LINE)

print(
    f"{'Stock Price':<35}: {stock_price:12.6f}"
)

print(
    f"{'Strike Price':<35}: {strike_price:12.6f}"
)

print(
    f"{'Interest Rate':<35}: {interest_rate * 100.0:11.6f}%"
)

print(
    f"{'Dividend Yield':<35}: {dividend_yield * 100.0:11.6f}%"
)

print(
    f"{'Volatility':<35}: {volatility * 100.0:11.6f}%"
)

print(
    f"{'CRR Steps':<35}: {num_steps:12d}"
)

eur_put_value = option_value(
    european_put,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    crr_model,
)

amer_put_value = option_value(
    american_put,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    crr_model,
)

eur_call_value = option_value(
    european_call,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    crr_model,
)

amer_call_value = option_value(
    american_call,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    crr_model,
)

print("\n" + SUBLINE)

print(
    f"{'OPTION':<30}"
    f"{'VALUE':>18}"
)

print(SUBLINE)

print(
    f"{'European Put':<30}"
    f"{eur_put_value:18.8f}"
)

print(
    f"{'American Put':<30}"
    f"{amer_put_value:18.8f}"
)

print(
    f"{'European Call':<30}"
    f"{eur_call_value:18.8f}"
)

print(
    f"{'American Call':<30}"
    f"{amer_call_value:18.8f}"
)


# ============================================================================
# 2. EARLY-EXERCISE PREMIUM
# ============================================================================
#
# An American option can be exercised before expiry.
#
# Therefore:
#
#       American Value >= European Value
#
# for otherwise identical contracts.
#
# The difference
#
#       American Value - European Value
#
# is the early-exercise premium.
# ============================================================================

print("\n" + LINE)
print("2. EARLY-EXERCISE PREMIUM")
print(LINE)

put_exercise_premium = (
    amer_put_value
    - eur_put_value
)

call_exercise_premium = (
    amer_call_value
    - eur_call_value
)

print(
    f"{'OPTION':<30}"
    f"{'EUROPEAN':>18}"
    f"{'AMERICAN':>18}"
    f"{'EXERCISE PREMIUM':>22}"
)

print(SUBLINE)

print(
    f"{'Put':<30}"
    f"{eur_put_value:18.8f}"
    f"{amer_put_value:18.8f}"
    f"{put_exercise_premium:22.8f}"
)

print(
    f"{'Call':<30}"
    f"{eur_call_value:18.8f}"
    f"{amer_call_value:18.8f}"
    f"{call_exercise_premium:22.8f}"
)


# ============================================================================
# 3. GREEKS
# ============================================================================

print("\n" + LINE)
print("3. OPTION GREEKS")
print(LINE)

print(
    f"{'OPTION':<25}"
    f"{'VALUE':>16}"
    f"{'DELTA':>16}"
    f"{'GAMMA':>16}"
    f"{'THETA':>16}"
)

print(SUBLINE)

for label, option in [
    ("European Put", european_put),
    ("American Put", american_put),
    ("European Call", european_call),
    ("American Call", american_call),
]:

    v = option.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        crr_model,
    )

    delta = option.delta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        crr_model,
    )

    gamma = option.gamma(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        crr_model,
    )

    theta = option.theta(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        crr_model,
    )

    print(
        f"{label:<25}"
        f"{v:16.8f}"
        f"{delta:16.8f}"
        f"{gamma:16.8f}"
        f"{theta:16.8f}"
    )


# ============================================================================
# 4. AVAILABLE FINANCEPY BLACK-SCHOLES MODEL TYPES
# ============================================================================

print("\n" + LINE)
print("4. AVAILABLE FINANCEPY BLACK-SCHOLES MODEL TYPES")
print(LINE)

for model_type in BlackScholesTypes:
    print(model_type)


# ============================================================================
# 5. AMERICAN OPTION PRICING METHOD COMPARISON
# ============================================================================
#
# A high-step CRR tree is used as the numerical benchmark.
#
# FinancePy analytic/semi-analytic model types that support American options
# are discovered dynamically and compared with the tree.
# ============================================================================

print("\n" + LINE)
print("5. AMERICAN OPTION PRICING METHOD COMPARISON")
print(LINE)

benchmark_steps = 2000

benchmark_model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    benchmark_steps,
)

benchmark_put = american_put.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

benchmark_call = american_call.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    benchmark_model,
)

print(
    f"{'CRR Benchmark Steps':<40}: "
    f"{benchmark_steps}"
)

print(
    f"{'CRR American Put':<40}: "
    f"{benchmark_put:12.8f}"
)

print(
    f"{'CRR American Call':<40}: "
    f"{benchmark_call:12.8f}"
)


# ============================================================================
# 5.1 ANALYTIC / SEMI-ANALYTIC APPROXIMATIONS
# ============================================================================

print("\n" + SUBLINE)
print("ANALYTIC / SEMI-ANALYTIC METHODS")
print(SUBLINE)

put_methods = available_american_models(
    american_put,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    volatility,
)

call_methods = available_american_models(
    american_call,
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    volatility,
)

print("\nAMERICAN PUT")

print(
    f"{'METHOD':<40}"
    f"{'VALUE':>16}"
    f"{'ERROR VS CRR':>18}"
)

print(SUBLINE)

for model_type, method_value in put_methods:

    error = (
        method_value
        - benchmark_put
    )

    print(
        f"{model_type_name(model_type):<40}"
        f"{method_value:16.8f}"
        f"{error:18.8f}"
    )

print("\nAMERICAN CALL")

print(
    f"{'METHOD':<40}"
    f"{'VALUE':>16}"
    f"{'ERROR VS CRR':>18}"
)

print(SUBLINE)

for model_type, method_value in call_methods:

    error = (
        method_value
        - benchmark_call
    )

    print(
        f"{model_type_name(model_type):<40}"
        f"{method_value:16.8f}"
        f"{error:18.8f}"
    )


# ============================================================================
# 6. CRR VERSUS LSMC
# ============================================================================

print("\n" + LINE)
print("6. CRR VERSUS LSMC")
print(LINE)

num_paths = 50000
lsmc_steps = 100

lsmc_model = BlackScholes(
    volatility,
    BlackScholesTypes.LSMC,
    lsmc_steps,
    num_paths,
)

start = time.time()

lsmc_put_value = american_put.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    lsmc_model,
)

lsmc_put_time = time.time() - start

start = time.time()

lsmc_call_value = american_call.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    lsmc_model,
)

lsmc_call_time = time.time() - start

print(
    f"{'METHOD':<25}"
    f"{'PUT VALUE':>18}"
    f"{'PUT ERROR':>18}"
    f"{'CALL VALUE':>18}"
    f"{'CALL ERROR':>18}"
)

print(SUBLINE)

print(
    f"{'CRR Benchmark':<25}"
    f"{benchmark_put:18.8f}"
    f"{0.0:18.8f}"
    f"{benchmark_call:18.8f}"
    f"{0.0:18.8f}"
)

print(
    f"{'LSMC':<25}"
    f"{lsmc_put_value:18.8f}"
    f"{lsmc_put_value - benchmark_put:18.8f}"
    f"{lsmc_call_value:18.8f}"
    f"{lsmc_call_value - benchmark_call:18.8f}"
)


# ============================================================================
# 7. OPTION VALUE VERSUS STOCK PRICE
# ============================================================================

print("\n" + LINE)
print("7. OPTION VALUE VERSUS STOCK PRICE")
print(LINE)

stock_prices = np.linspace(
    20.0,
    80.0,
    61,
)

european_put_values = []
american_put_values = []

european_call_values = []
american_call_values = []

put_intrinsic_values = []
call_intrinsic_values = []

plot_model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    500,
)

for s in stock_prices:

    european_put_values.append(
        european_put.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    american_put_values.append(
        american_put.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    european_call_values.append(
        european_call.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    american_call_values.append(
        american_call.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    put_intrinsic_values.append(
        max(
            strike_price - s,
            0.0,
        )
    )

    call_intrinsic_values.append(
        max(
            s - strike_price,
            0.0,
        )
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
    label="Put Intrinsic Value",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("European and American Put Value")
plt.grid(True)
plt.legend()
plt.show()


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
    label="Call Intrinsic Value",
)

plt.xlabel("Stock Price")
plt.ylabel("Option Value")
plt.title("European and American Call Value")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. EARLY-EXERCISE PREMIUM VERSUS STOCK PRICE
# ============================================================================

print("\n" + LINE)
print("8. EARLY-EXERCISE PREMIUM VERSUS STOCK PRICE")
print(LINE)

put_exercise_premiums = (
    np.asarray(american_put_values)
    - np.asarray(european_put_values)
)

call_exercise_premiums = (
    np.asarray(american_call_values)
    - np.asarray(european_call_values)
)

plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    put_exercise_premiums,
    label="American Put Premium",
)

plt.plot(
    stock_prices,
    call_exercise_premiums,
    label="American Call Premium",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Stock Price")
plt.ylabel("American - European Value")
plt.title("Early-Exercise Premium")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. ANALYTIC APPROXIMATION ERROR VERSUS STOCK PRICE
# ============================================================================
#
# The CRR tree is used as the benchmark.
#
# Plotting approximation error is generally more informative than plotting
# several almost identical option-value curves.
# ============================================================================

print("\n" + LINE)
print("9. ANALYTIC APPROXIMATION ERROR VERSUS STOCK PRICE")
print(LINE)

for model_type, _ in put_methods:

    errors = []

    approximation_model = BlackScholes(
        volatility,
        model_type,
    )

    for s in stock_prices:

        approximation_value = american_put.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            approximation_model,
        )

        benchmark_value = american_put.value(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )

        errors.append(
            approximation_value
            - benchmark_value
        )

    plt.figure(figsize=(9, 6))

    plt.plot(
        stock_prices,
        errors,
    )

    plt.axhline(
        0.0,
        linestyle="--",
    )

    plt.xlabel("Stock Price")
    plt.ylabel("Approximation - CRR")
    plt.title(
        "American Put Approximation Error: "
        + model_type_name(model_type)
    )

    plt.grid(True)
    plt.show()


# ============================================================================
# 10. EFFECT OF DIVIDENDS ON THE AMERICAN CALL
# ============================================================================
#
# For a non-dividend-paying stock, early exercise of a standard American call
# is not optimal under the usual Black-Scholes assumptions.
#
# Dividends change this trade-off because exercising the option converts the
# option into stock ownership and therefore gives access to future dividends.
# ============================================================================

print("\n" + LINE)
print("10. EFFECT OF DIVIDENDS ON THE AMERICAN CALL")
print(LINE)

dividend_yields = np.linspace(
    0.0,
    0.12,
    25,
)

european_values = []
american_values = []
exercise_premiums = []

for q in dividend_yields:

    _, q_curve = build_curves(
        value_dt,
        interest_rate,
        q,
    )

    european_value = european_call.value(
        value_dt,
        stock_price,
        discount_curve,
        q_curve,
        plot_model,
    )

    american_value = american_call.value(
        value_dt,
        stock_price,
        discount_curve,
        q_curve,
        plot_model,
    )

    european_values.append(
        european_value
    )

    american_values.append(
        american_value
    )

    exercise_premiums.append(
        american_value
        - european_value
    )


plt.figure(figsize=(9, 6))

plt.plot(
    dividend_yields * 100.0,
    european_values,
    label="European Call",
)

plt.plot(
    dividend_yields * 100.0,
    american_values,
    label="American Call",
)

plt.xlabel("Dividend Yield (%)")
plt.ylabel("Option Value")
plt.title("Call Value Versus Dividend Yield")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(9, 6))

plt.plot(
    dividend_yields * 100.0,
    exercise_premiums,
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Dividend Yield (%)")
plt.ylabel("American - European Value")
plt.title("American Call Early-Exercise Premium")
plt.grid(True)
plt.show()


# ============================================================================
# 11. OPTION VALUE VERSUS VOLATILITY
# ============================================================================

print("\n" + LINE)
print("11. OPTION VALUE VERSUS VOLATILITY")
print(LINE)

volatilities = np.linspace(
    0.10,
    0.80,
    36,
)

put_vol_values = []
call_vol_values = []

for vol in volatilities:

    vol_model = BlackScholes(
        vol,
        BlackScholesTypes.CRR_TREE,
        500,
    )

    put_vol_values.append(
        american_put.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            vol_model,
        )
    )

    call_vol_values.append(
        american_call.value(
            value_dt,
            stock_price,
            discount_curve,
            dividend_curve,
            vol_model,
        )
    )


plt.figure(figsize=(9, 6))

plt.plot(
    volatilities * 100.0,
    put_vol_values,
    label="American Put",
)

plt.plot(
    volatilities * 100.0,
    call_vol_values,
    label="American Call",
)

plt.xlabel("Volatility (%)")
plt.ylabel("Option Value")
plt.title("American Option Value Versus Volatility")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 12. OPTION VALUE AS EXPIRY APPROACHES
# ============================================================================
#
# Keep the final expiry date fixed and move the valuation date forward.
# This illustrates the loss of time value as the option approaches expiry.
# ============================================================================

print("\n" + LINE)
print("12. OPTION VALUE AS EXPIRY APPROACHES")
print(LINE)

days_to_expiry = np.array(
    [
        365,
        300,
        240,
        180,
        120,
        90,
        60,
        30,
        14,
        7,
        2,
        1,
    ]
)

time_put_values = []
time_call_values = []

for days in days_to_expiry:

    moving_value_dt = expiry_dt.add_days(
        -int(days)
    )

    moving_discount_curve = FlatDiscountCurve(
        moving_value_dt,
        interest_rate,
    )

    moving_dividend_curve = FlatDiscountCurve(
        moving_value_dt,
        dividend_yield,
    )

    moving_model = BlackScholes(
        volatility,
        BlackScholesTypes.CRR_TREE,
        500,
    )

    time_put_values.append(
        american_put.value(
            moving_value_dt,
            stock_price,
            moving_discount_curve,
            moving_dividend_curve,
            moving_model,
        )
    )

    time_call_values.append(
        american_call.value(
            moving_value_dt,
            stock_price,
            moving_discount_curve,
            moving_dividend_curve,
            moving_model,
        )
    )


plt.figure(figsize=(9, 6))

plt.plot(
    days_to_expiry,
    time_put_values,
    marker="o",
    label="American Put",
)

plt.plot(
    days_to_expiry,
    time_call_values,
    marker="o",
    label="American Call",
)

plt.xlabel("Days to Expiry")
plt.ylabel("Option Value")
plt.title("American Option Value as Expiry Approaches")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 13. DELTA VERSUS STOCK PRICE
# ============================================================================

print("\n" + LINE)
print("13. DELTA VERSUS STOCK PRICE")
print(LINE)

put_deltas = []
call_deltas = []

for s in stock_prices:

    put_deltas.append(
        american_put.delta(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    call_deltas.append(
        american_call.delta(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )


plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    put_deltas,
    label="American Put Delta",
)

plt.plot(
    stock_prices,
    call_deltas,
    label="American Call Delta",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Stock Price")
plt.ylabel("Delta")
plt.title("American Option Delta")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 14. GAMMA VERSUS STOCK PRICE
# ============================================================================

print("\n" + LINE)
print("14. GAMMA VERSUS STOCK PRICE")
print(LINE)

put_gammas = []
call_gammas = []

for s in stock_prices:

    put_gammas.append(
        american_put.gamma(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )

    call_gammas.append(
        american_call.gamma(
            value_dt,
            s,
            discount_curve,
            dividend_curve,
            plot_model,
        )
    )


plt.figure(figsize=(9, 6))

plt.plot(
    stock_prices,
    put_gammas,
    label="American Put Gamma",
)

plt.plot(
    stock_prices,
    call_gammas,
    label="American Call Gamma",
)

plt.xlabel("Stock Price")
plt.ylabel("Gamma")
plt.title("American Option Gamma")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 15. CRR TREE CONVERGENCE
# ============================================================================
#
# A binomial tree converges toward the continuous-time American-option value
# as the number of time steps increases. Binomial convergence is generally
# not perfectly monotonic because odd/even tree structures can oscillate.
# ============================================================================

print("\n" + LINE)
print("15. CRR TREE CONVERGENCE")
print(LINE)

tree_steps = np.array(
    [
        25,
        50,
        100,
        200,
        500,
        1000,
        2000,
    ]
)

tree_put_values = []
tree_call_values = []
tree_times = []

print(
    f"{'STEPS':>10}"
    f"{'PUT':>18}"
    f"{'CALL':>18}"
    f"{'TIME':>18}"
)

print(SUBLINE)

for steps in tree_steps:

    convergence_model = BlackScholes(
        volatility,
        BlackScholesTypes.CRR_TREE,
        int(steps),
    )

    start = time.time()

    put_value = american_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        convergence_model,
    )

    call_value = american_call.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        convergence_model,
    )

    elapsed = time.time() - start

    tree_put_values.append(
        put_value
    )

    tree_call_values.append(
        call_value
    )

    tree_times.append(
        elapsed
    )

    print(
        f"{steps:10d}"
        f"{put_value:18.8f}"
        f"{call_value:18.8f}"
        f"{elapsed:18.8f}"
    )


plt.figure(figsize=(9, 6))

plt.plot(
    tree_steps,
    tree_put_values,
    marker="o",
    label="American Put",
)

plt.axhline(
    benchmark_put,
    linestyle="--",
    label="2000-Step Benchmark",
)

plt.xlabel("Number of CRR Steps")
plt.ylabel("Option Value")
plt.title("CRR Tree Convergence - American Put")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 16. COMPUTATIONAL SPEED COMPARISON
# ============================================================================

print("\n" + LINE)
print("16. COMPUTATIONAL SPEED COMPARISON")
print(LINE)

print(
    f"{'METHOD':<40}"
    f"{'VALUE':>18}"
    f"{'TIME':>18}"
)

print(SUBLINE)

for steps in [
    100,
    500,
    1000,
    2000,
]:

    speed_model = BlackScholes(
        volatility,
        BlackScholesTypes.CRR_TREE,
        steps,
    )

    start = time.time()

    method_value = american_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        speed_model,
    )

    elapsed = time.time() - start

    print(
        f"{'CRR ' + str(steps):<40}"
        f"{method_value:18.8f}"
        f"{elapsed:18.8f}"
    )


for model_type, _ in put_methods:

    speed_model = BlackScholes(
        volatility,
        model_type,
    )

    start = time.time()

    method_value = american_put.value(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        speed_model,
    )

    elapsed = time.time() - start

    print(
        f"{model_type_name(model_type):<40}"
        f"{method_value:18.8f}"
        f"{elapsed:18.8f}"
    )


start = time.time()

method_value = american_put.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    lsmc_model,
)

elapsed = time.time() - start

print(
    f"{'LSMC':<40}"
    f"{method_value:18.8f}"
    f"{elapsed:18.8f}"
)


# ============================================================================
# 17. DELTA BUMP-AND-REVALUE TEST
# ============================================================================
#
# Delta is:
#
#                   dV
#           Delta = --
#                   dS
#
# A central finite-difference approximation is:
#
#              V(S + dS) - V(S - dS)
#     Delta ~= -----------------------
#                       2 dS
#
# This provides an independent numerical check of option.delta().
# ============================================================================

print("\n" + LINE)
print("17. DELTA BUMP-AND-REVALUE TEST")
print(LINE)

test_option = american_put

test_model = BlackScholes(
    volatility,
    BlackScholesTypes.CRR_TREE,
    2000,
)

stock_bump = 0.01

base_value = test_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    test_model,
)

value_up = test_option.value(
    value_dt,
    stock_price + stock_bump,
    discount_curve,
    dividend_curve,
    test_model,
)

value_down = test_option.value(
    value_dt,
    stock_price - stock_bump,
    discount_curve,
    dividend_curve,
    test_model,
)

delta_bump = (
    value_up
    - value_down
) / (
    2.0 * stock_bump
)

delta_function = test_option.delta(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    test_model,
)

delta_difference = (
    delta_function
    - delta_bump
)

print(
    f"{'Base Value':<40}: "
    f"{base_value:16.10f}"
)

print(
    f"{'Stock Bump':<40}: "
    f"{stock_bump:16.10f}"
)

print(
    f"{'Value S + dS':<40}: "
    f"{value_up:16.10f}"
)

print(
    f"{'Value S - dS':<40}: "
    f"{value_down:16.10f}"
)

print(
    f"{'FinancePy Delta':<40}: "
    f"{delta_function:16.10f}"
)

print(
    f"{'Bump-and-Revalue Delta':<40}: "
    f"{delta_bump:16.10f}"
)

print(
    f"{'Difference':<40}: "
    f"{delta_difference:16.10f}"
)


# ============================================================================
# 18. GAMMA BUMP-AND-REVALUE TEST
# ============================================================================
#
# Gamma measures the rate of change of delta with stock price:
#
#                    d2V
#           Gamma = -----
#                    dS2
#
# The central finite-difference approximation is:
#
#          V(S+dS) - 2V(S) + V(S-dS)
# Gamma ~= ---------------------------
#                       dS^2
# ============================================================================

print("\n" + LINE)
print("18. GAMMA BUMP-AND-REVALUE TEST")
print(LINE)

gamma_bump = (
    value_up
    - 2.0 * base_value
    + value_down
) / (
    stock_bump * stock_bump
)

gamma_function = test_option.gamma(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    test_model,
)

gamma_difference = (
    gamma_function
    - gamma_bump
)

print(
    f"{'FinancePy Gamma':<40}: "
    f"{gamma_function:16.10f}"
)

print(
    f"{'Bump-and-Revalue Gamma':<40}: "
    f"{gamma_bump:16.10f}"
)

print(
    f"{'Difference':<40}: "
    f"{gamma_difference:16.10f}"
)


# ============================================================================
# 19. SUMMARY
# ============================================================================

print("\n" + LINE)
print("19. SUMMARY")
print(LINE)

print(
    """
This example illustrates several important properties of American options.

1. EUROPEAN VERSUS AMERICAN

   An American option has all of the exercise opportunities of the equivalent
   European option plus the ability to exercise before expiry. Its value
   should therefore not be below the corresponding European option value.

2. EARLY-EXERCISE PREMIUM

   American Value - European Value

   measures the value of the additional early-exercise right.

3. AMERICAN PUTS

   Early exercise can be valuable for puts. A sufficiently deep in-the-money
   put may be worth exercising early because the strike proceeds can then be
   received and invested immediately.

4. AMERICAN CALLS

   For a non-dividend-paying stock, early exercise of a standard call is not
   optimal under the usual Black-Scholes assumptions.

   Positive dividends can make early exercise valuable because exercising the
   option gives ownership of the stock and therefore access to dividends.

5. CRR TREE

   The CRR tree handles early exercise explicitly by comparing continuation
   value with exercise value at every node.

   Increasing the number of steps generally improves accuracy, although
   binomial prices can oscillate as the number of steps changes.

6. LSMC

   Longstaff-Schwartz Monte Carlo estimates continuation values using
   regression. It is particularly useful when the state space becomes too
   large for simple lattice methods.

7. ANALYTIC / SEMI-ANALYTIC APPROXIMATIONS

   These methods seek to approximate the American early-exercise feature much
   faster than a large tree or Monte Carlo simulation.

   Their errors should therefore be studied over a range of stock prices and
   market conditions rather than at only one point.

8. DELTA

   Delta measures first-order sensitivity to the stock price.

   The bump-and-revalue calculation provides an independent numerical check
   of FinancePy's delta calculation.

9. GAMMA

   Gamma measures the curvature of option value with respect to stock price.

   The central second difference provides an independent numerical check of
   FinancePy's gamma calculation.

10. MODEL CHOICE

   A useful practical comparison is therefore:

       analytic/semi-analytic approximation
               versus
       CRR tree
               versus
       LSMC

   in terms of both valuation accuracy and computational cost.
"""
)

print(LINE)


# ============================================================================
# OPTIONAL: REPRODUCE LONGSTAFF-SCHWARTZ PAPER EXAMPLE
# ============================================================================

# replicate_ls_paper()
