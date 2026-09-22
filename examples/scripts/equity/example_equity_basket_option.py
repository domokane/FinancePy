# ============================================================================
# FINANCEPY EXAMPLES - EquityBasketOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_basket_option import EquityBasketOption
from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style
from financepy.utils.global_types import OptionTypes
from financepy.utils.helpers import beta_vector_to_corr_matrix


set_plot_style()

LINE = "=" * 78


# ============================================================================
# SUPPORTING FUNCTIONS
# ============================================================================


def build_dividend_curves(value_dt, dividend_yields):
    """Build one flat dividend curve for each asset."""

    dividend_curves = []

    for dividend_yield in dividend_yields:

        dividend_curve = FlatDiscountCurve(
            value_dt,
            dividend_yield,
        )

        dividend_curves.append(
            dividend_curve,
        )

    return dividend_curves


def build_correlation_matrix(num_assets, beta):
    """Build the one-factor correlation matrix from equal asset betas."""

    betas = np.ones(num_assets) * beta

    corr_matrix = beta_vector_to_corr_matrix(
        betas,
    )

    return corr_matrix


def pairwise_correlation(corr_matrix):
    """Return the common off-diagonal correlation for a homogeneous beta."""

    return corr_matrix[0, 1]


def value_basket(
    option,
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
    num_paths,
):
    """Return analytic and Monte Carlo basket-option values."""

    analytic_value = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    mc_value = option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    return analytic_value, mc_value


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

interest_rate = 0.05

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

num_assets = 5
strike_price = 100.0

num_paths = 10000

beta_list = np.linspace(
    0.0,
    0.999999,
    11,
)


# ============================================================================
# 1. HOMOGENEOUS BASKET CALL
# ============================================================================
# What this section demonstrates:
#
# Values a European call on a homogeneous equity basket.
#
# Every asset has the same:
#
#   - stock price
#   - volatility
#   - dividend yield
#
# The common factor loading beta is varied. The corresponding pairwise
# asset correlation is obtained from the resulting correlation matrix.
#
# The analytic basket approximation is compared with Monte Carlo.
# ============================================================================

print("\n" + LINE)
print("1. HOMOGENEOUS BASKET CALL")
print(LINE)

stock_prices = np.ones(num_assets) * 100.0
volatilities = np.ones(num_assets) * 0.30
dividend_yields = np.ones(num_assets) * 0.01

dividend_curves = build_dividend_curves(
    value_dt,
    dividend_yields,
)

call_option = EquityBasketOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    num_assets,
)

hom_call_correlations = []
hom_call_analytic_values = []
hom_call_mc_values = []

print(
    f"{'BETA':>10}"
    f"{'CORRELATION':>15}"
    f"{'ANALYTIC':>15}"
    f"{'MONTE CARLO':>15}"
    f"{'DIFFERENCE':>15}"
    f"{'TIME':>12}"
)

print("-" * 82)

for beta in beta_list:

    corr_matrix = build_correlation_matrix(
        num_assets,
        beta,
    )

    correlation = pairwise_correlation(
        corr_matrix,
    )

    start = time.time()

    analytic_value, mc_value = value_basket(
        call_option,
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    end = time.time()

    duration = end - start
    difference = mc_value - analytic_value

    hom_call_correlations.append(
        correlation,
    )

    hom_call_analytic_values.append(
        analytic_value,
    )

    hom_call_mc_values.append(
        mc_value,
    )

    print(
        f"{beta:10.4f}"
        f"{correlation:15.4f}"
        f"{analytic_value:15.6f}"
        f"{mc_value:15.6f}"
        f"{difference:15.6f}"
        f"{duration:12.4f}"
    )

hom_call_correlations = np.array(
    hom_call_correlations,
)

hom_call_analytic_values = np.array(
    hom_call_analytic_values,
)

hom_call_mc_values = np.array(
    hom_call_mc_values,
)


# ============================================================================
# 2. PLOT HOMOGENEOUS CALL VALUE VERSUS CORRELATION
# ============================================================================
# What this section demonstrates:
#
# Correlation is an important pricing input for a basket option.
#
# When correlation is low, movements in the individual stocks diversify one
# another. As correlation increases, the assets behave increasingly like a
# common underlying.
#
# The graph also compares the analytic basket approximation with Monte Carlo.
# ============================================================================

print("\n" + LINE)
print("2. HOMOGENEOUS CALL VALUE VERSUS CORRELATION")
print(LINE)

plt.figure()

plt.plot(
    hom_call_correlations,
    hom_call_analytic_values,
    marker="o",
    label="Analytic",
)

plt.plot(
    hom_call_correlations,
    hom_call_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Pairwise Asset Correlation")
plt.ylabel("Basket Call Value")

plt.title(
    "Homogeneous Basket Call Value versus Correlation"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. HETEROGENEOUS BASKET CALL
# ============================================================================
# What this section demonstrates:
#
# Repeats the call-option experiment for a heterogeneous basket.
#
# The constituent assets now have different:
#
#   - stock prices
#   - volatilities
#   - dividend yields
#
# This provides a more realistic test of the basket approximation.
# ============================================================================

print("\n" + LINE)
print("3. HETEROGENEOUS BASKET CALL")
print(LINE)

stock_prices = np.array(
    [
        100.0,
        105.0,
        120.0,
        100.0,
        90.0,
    ]
)

volatilities = np.array(
    [
        0.30,
        0.20,
        0.25,
        0.22,
        0.40,
    ]
)

dividend_yields = np.array(
    [
        0.01,
        0.02,
        0.04,
        0.01,
        0.02,
    ]
)

dividend_curves = build_dividend_curves(
    value_dt,
    dividend_yields,
)

hetero_call_correlations = []
hetero_call_analytic_values = []
hetero_call_mc_values = []

print(
    f"{'BETA':>10}"
    f"{'CORRELATION':>15}"
    f"{'ANALYTIC':>15}"
    f"{'MONTE CARLO':>15}"
    f"{'DIFFERENCE':>15}"
)

print("-" * 70)

for beta in beta_list:

    corr_matrix = build_correlation_matrix(
        num_assets,
        beta,
    )

    correlation = pairwise_correlation(
        corr_matrix,
    )

    analytic_value, mc_value = value_basket(
        call_option,
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    difference = mc_value - analytic_value

    hetero_call_correlations.append(
        correlation,
    )

    hetero_call_analytic_values.append(
        analytic_value,
    )

    hetero_call_mc_values.append(
        mc_value,
    )

    print(
        f"{beta:10.4f}"
        f"{correlation:15.4f}"
        f"{analytic_value:15.6f}"
        f"{mc_value:15.6f}"
        f"{difference:15.6f}"
    )

hetero_call_correlations = np.array(
    hetero_call_correlations,
)

hetero_call_analytic_values = np.array(
    hetero_call_analytic_values,
)

hetero_call_mc_values = np.array(
    hetero_call_mc_values,
)


# ============================================================================
# 4. HOMOGENEOUS BASKET PUT
# ============================================================================
# What this section demonstrates:
#
# Repeats the homogeneous-basket experiment for a European put.
#
# This checks that the effect of correlation and the agreement between the
# analytic approximation and Monte Carlo are not specific to call options.
# ============================================================================

print("\n" + LINE)
print("4. HOMOGENEOUS BASKET PUT")
print(LINE)

stock_prices = np.ones(num_assets) * 100.0
volatilities = np.ones(num_assets) * 0.30
dividend_yields = np.ones(num_assets) * 0.01

dividend_curves = build_dividend_curves(
    value_dt,
    dividend_yields,
)

put_option = EquityBasketOption(
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_PUT,
    num_assets,
)

hom_put_correlations = []
hom_put_analytic_values = []
hom_put_mc_values = []

print(
    f"{'BETA':>10}"
    f"{'CORRELATION':>15}"
    f"{'ANALYTIC':>15}"
    f"{'MONTE CARLO':>15}"
    f"{'DIFFERENCE':>15}"
)

print("-" * 70)

for beta in beta_list:

    corr_matrix = build_correlation_matrix(
        num_assets,
        beta,
    )

    correlation = pairwise_correlation(
        corr_matrix,
    )

    analytic_value, mc_value = value_basket(
        put_option,
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    difference = mc_value - analytic_value

    hom_put_correlations.append(
        correlation,
    )

    hom_put_analytic_values.append(
        analytic_value,
    )

    hom_put_mc_values.append(
        mc_value,
    )

    print(
        f"{beta:10.4f}"
        f"{correlation:15.4f}"
        f"{analytic_value:15.6f}"
        f"{mc_value:15.6f}"
        f"{difference:15.6f}"
    )

hom_put_correlations = np.array(
    hom_put_correlations,
)

hom_put_analytic_values = np.array(
    hom_put_analytic_values,
)

hom_put_mc_values = np.array(
    hom_put_mc_values,
)


# ============================================================================
# 5. HETEROGENEOUS BASKET PUT
# ============================================================================
# What this section demonstrates:
#
# Values a European put on the heterogeneous basket used previously for the
# call option.
#
# This completes the analytic-versus-Monte-Carlo comparison for both calls
# and puts under homogeneous and heterogeneous market inputs.
# ============================================================================

print("\n" + LINE)
print("5. HETEROGENEOUS BASKET PUT")
print(LINE)

stock_prices = np.array(
    [
        100.0,
        105.0,
        120.0,
        100.0,
        90.0,
    ]
)

volatilities = np.array(
    [
        0.30,
        0.20,
        0.25,
        0.22,
        0.40,
    ]
)

dividend_yields = np.array(
    [
        0.01,
        0.02,
        0.04,
        0.01,
        0.02,
    ]
)

dividend_curves = build_dividend_curves(
    value_dt,
    dividend_yields,
)

hetero_put_correlations = []
hetero_put_analytic_values = []
hetero_put_mc_values = []

print(
    f"{'BETA':>10}"
    f"{'CORRELATION':>15}"
    f"{'ANALYTIC':>15}"
    f"{'MONTE CARLO':>15}"
    f"{'DIFFERENCE':>15}"
)

print("-" * 70)

for beta in beta_list:

    corr_matrix = build_correlation_matrix(
        num_assets,
        beta,
    )

    correlation = pairwise_correlation(
        corr_matrix,
    )

    analytic_value, mc_value = value_basket(
        put_option,
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
    )

    difference = mc_value - analytic_value

    hetero_put_correlations.append(
        correlation,
    )

    hetero_put_analytic_values.append(
        analytic_value,
    )

    hetero_put_mc_values.append(
        mc_value,
    )

    print(
        f"{beta:10.4f}"
        f"{correlation:15.4f}"
        f"{analytic_value:15.6f}"
        f"{mc_value:15.6f}"
        f"{difference:15.6f}"
    )

hetero_put_correlations = np.array(
    hetero_put_correlations,
)

hetero_put_analytic_values = np.array(
    hetero_put_analytic_values,
)

hetero_put_mc_values = np.array(
    hetero_put_mc_values,
)


# ============================================================================
# 6. ANALYTIC APPROXIMATION ERROR
# ============================================================================
# What this section demonstrates:
#
# Plots:
#
#                   Monte Carlo - Analytic
#
# against correlation.
#
# Monte Carlo provides an independent numerical benchmark for assessing the
# basket-option analytic approximation.
#
# The error need not be exactly zero because:
#
#   - the analytic method is an approximation
#   - Monte Carlo contains sampling error
#
# Comparing homogeneous and heterogeneous baskets shows whether approximation
# quality changes when the constituent assets have different characteristics.
# ============================================================================

print("\n" + LINE)
print("6. ANALYTIC APPROXIMATION ERROR")
print(LINE)

hom_call_errors = (
    hom_call_mc_values
    - hom_call_analytic_values
)

hetero_call_errors = (
    hetero_call_mc_values
    - hetero_call_analytic_values
)

plt.figure()

plt.plot(
    hom_call_correlations,
    hom_call_errors,
    marker="o",
    label="Homogeneous Call",
)

plt.plot(
    hetero_call_correlations,
    hetero_call_errors,
    marker="o",
    label="Heterogeneous Call",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Pairwise Asset Correlation")
plt.ylabel("Monte Carlo - Analytic")

plt.title(
    "Basket Option Analytic Approximation Error"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. CORRELATION AND DIVERSIFICATION
# ============================================================================
# What this section demonstrates:
#
# Compares homogeneous call and put values as correlation changes.
#
# Correlation controls the amount of diversification within the basket.
#
# Low correlation allows individual stock movements to offset one another,
# reducing the volatility of the basket.
#
# High correlation reduces this diversification benefit because the assets
# increasingly move together.
#
# Basket option values therefore depend materially on correlation even though
# the individual stock volatilities remain unchanged.
# ============================================================================

print("\n" + LINE)
print("7. CORRELATION AND DIVERSIFICATION")
print(LINE)

plt.figure()

plt.plot(
    hom_call_correlations,
    hom_call_analytic_values,
    marker="o",
    label="Call",
)

plt.plot(
    hom_put_correlations,
    hom_put_analytic_values,
    marker="o",
    label="Put",
)

plt.xlabel("Pairwise Asset Correlation")
plt.ylabel("Basket Option Value")

plt.title(
    "Basket Option Value and Diversification"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. MONTE CARLO CONVERGENCE
# ============================================================================
# What this section demonstrates:
#
# Holds the basket and correlation fixed while increasing the number of
# Monte Carlo paths.
#
# As the number of paths increases, sampling noise should generally decrease.
#
# The analytic approximation is shown as a reference value. The Monte Carlo
# calculation is not expected to converge exactly to the analytic value if
# the analytic method itself contains approximation error.
# ============================================================================

print("\n" + LINE)
print("8. MONTE CARLO CONVERGENCE")
print(LINE)

stock_prices = np.ones(num_assets) * 100.0
volatilities = np.ones(num_assets) * 0.30
dividend_yields = np.ones(num_assets) * 0.01

dividend_curves = build_dividend_curves(
    value_dt,
    dividend_yields,
)

beta = np.sqrt(0.50)

corr_matrix = build_correlation_matrix(
    num_assets,
    beta,
)

correlation = pairwise_correlation(
    corr_matrix,
)

analytic_value = call_option.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
)

path_counts = np.array(
    [
        1000,
        2000,
        5000,
        10000,
        20000,
        50000,
    ]
)

mc_convergence_values = []
mc_convergence_errors = []

print(
    f"{'NUM PATHS':>15}"
    f"{'MC VALUE':>15}"
    f"{'ANALYTIC':>15}"
    f"{'DIFFERENCE':>15}"
    f"{'TIME':>15}"
)

print("-" * 75)

for path_count in path_counts:

    start = time.time()

    mc_value = call_option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        int(path_count),
    )

    end = time.time()

    duration = end - start

    error = (
        mc_value
        - analytic_value
    )

    mc_convergence_values.append(
        mc_value,
    )

    mc_convergence_errors.append(
        error,
    )

    print(
        f"{path_count:15d}"
        f"{mc_value:15.6f}"
        f"{analytic_value:15.6f}"
        f"{error:15.6f}"
        f"{duration:15.4f}"
    )

mc_convergence_values = np.array(
    mc_convergence_values,
)

mc_convergence_errors = np.array(
    mc_convergence_errors,
)

print(
    f"\nPairwise asset correlation = {correlation:.4f}"
)

plt.figure()

plt.plot(
    path_counts,
    mc_convergence_values,
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
plt.ylabel("Basket Call Value")

plt.title(
    "Basket Option Monte Carlo Convergence"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. MONTE CARLO ERROR VERSUS NUMBER OF PATHS
# ============================================================================
# What this section demonstrates:
#
# Shows the difference between the Monte Carlo value and the analytic
# approximation as the number of simulation paths increases.
#
# This graph makes the sampling behaviour easier to see than plotting the
# two option values directly.
#
# A persistent non-zero difference at large path counts may reflect the
# approximation error of the analytic basket model rather than Monte Carlo
# sampling error.
# ============================================================================

print("\n" + LINE)
print("9. MONTE CARLO ERROR VERSUS NUMBER OF PATHS")
print(LINE)

plt.figure()

plt.plot(
    path_counts,
    mc_convergence_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xscale("log")

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Monte Carlo - Analytic")

plt.title(
    "Basket Option Monte Carlo Error"
)

plt.grid(True)
plt.show()
