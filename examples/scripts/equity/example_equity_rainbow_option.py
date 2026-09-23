# ============================================================================
# FINANCEPY EXAMPLES - EquityRainbowOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

from math import sqrt
import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.global_types import OptionTypes
from financepy.utils.date import Date
from financepy.utils.helpers import beta_vector_to_corr_matrix

from financepy.products.equity.equity_vanilla_option import EquityVanillaOption
from financepy.models.black_scholes import BlackScholes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_rainbow_option import (
    EquityRainbowOption,
    EquityRainbowOptionTypes,
)


def print_header(number, title):
    """Print a section header."""

    print("\n" + "=" * 78)
    print(f"{number}. {title}")
    print("=" * 78)


def make_corr_matrix(num_assets, correlation):
    """Create an equicorrelation matrix."""

    betas = np.ones(num_assets) * sqrt(correlation)
    return beta_vector_to_corr_matrix(betas)


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

interest_rate = 0.05
dividend_yield = 0.01

num_assets = 2

stock_prices = np.ones(num_assets) * 100.0
volatilities = np.ones(num_assets) * 0.30

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curves = [
    FlatDiscountCurve(
        value_dt,
        dividend_yield,
    )
    for _ in range(num_assets)
]

strike = 100.0

num_paths = 100000
seed = 4242

corr_list = [
    0.00,
    0.10,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    0.99,
]


# ============================================================================
# 1. CALL ON MAXIMUM
# ============================================================================

print_header(1, "CALL ON MAXIMUM")

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

call_max_analytic = []
call_max_mc = []

print(
    f"{'CORRELATION':>12}"
    f"{'ANALYTIC':>16}"
    f"{'MC':>16}"
    f"{'MC-ANALYTIC':>16}"
)

print("-" * 60)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    value = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
        seed=seed,
    )

    call_max_analytic.append(value)
    call_max_mc.append(value_mc)

    print(
        f"{correlation:12.4f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{value_mc - value:16.8f}"
    )


# ============================================================================
# 2. CALL ON MINIMUM
# ============================================================================

print_header(2, "CALL ON MINIMUM")

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    [strike],
    num_assets,
)

call_min_analytic = []
call_min_mc = []

print(
    f"{'CORRELATION':>12}"
    f"{'ANALYTIC':>16}"
    f"{'MC':>16}"
    f"{'MC-ANALYTIC':>16}"
)

print("-" * 60)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    value = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
        seed=seed,
    )

    call_min_analytic.append(value)
    call_min_mc.append(value_mc)

    print(
        f"{correlation:12.4f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{value_mc - value:16.8f}"
    )


# ============================================================================
# 3. PUT ON MAXIMUM
# ============================================================================

print_header(3, "PUT ON MAXIMUM")

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MAXIMUM,
    [strike],
    num_assets,
)

put_max_analytic = []
put_max_mc = []

print(
    f"{'CORRELATION':>12}"
    f"{'ANALYTIC':>16}"
    f"{'MC':>16}"
    f"{'MC-ANALYTIC':>16}"
)

print("-" * 60)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    value = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
        seed=seed,
    )

    put_max_analytic.append(value)
    put_max_mc.append(value_mc)

    print(
        f"{correlation:12.4f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{value_mc - value:16.8f}"
    )


# ============================================================================
# 4. PUT ON MINIMUM
# ============================================================================

print_header(4, "PUT ON MINIMUM")

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MINIMUM,
    [strike],
    num_assets,
)

put_min_analytic = []
put_min_mc = []

print(
    f"{'CORRELATION':>12}"
    f"{'ANALYTIC':>16}"
    f"{'MC':>16}"
    f"{'MC-ANALYTIC':>16}"
)

print("-" * 60)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    value = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    value_mc = option.value_mc(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths,
        seed=seed,
    )

    put_min_analytic.append(value)
    put_min_mc.append(value_mc)

    print(
        f"{correlation:12.4f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{value_mc - value:16.8f}"
    )


# ============================================================================
# 5. OPTION VALUE VERSUS CORRELATION
# ============================================================================

print_header(5, "OPTION VALUE VERSUS CORRELATION")

fig, ax = plt.subplots()

ax.plot(
    corr_list,
    call_max_analytic,
    marker="o",
    label="Call max analytic",
)

ax.plot(
    corr_list,
    call_max_mc,
    marker="x",
    linestyle="--",
    label="Call max MC",
)

ax.plot(
    corr_list,
    call_min_analytic,
    marker="o",
    label="Call min analytic",
)

ax.plot(
    corr_list,
    call_min_mc,
    marker="x",
    linestyle="--",
    label="Call min MC",
)

ax.set_xlabel("Correlation")
ax.set_ylabel("Option value")
ax.set_title("Rainbow Call Value versus Correlation")
ax.grid(True)
ax.legend()

plt.show()


fig, ax = plt.subplots()

ax.plot(
    corr_list,
    put_max_analytic,
    marker="o",
    label="Put max analytic",
)

ax.plot(
    corr_list,
    put_max_mc,
    marker="x",
    linestyle="--",
    label="Put max MC",
)

ax.plot(
    corr_list,
    put_min_analytic,
    marker="o",
    label="Put min analytic",
)

ax.plot(
    corr_list,
    put_min_mc,
    marker="x",
    linestyle="--",
    label="Put min MC",
)

ax.set_xlabel("Correlation")
ax.set_ylabel("Option value")
ax.set_title("Rainbow Put Value versus Correlation")
ax.grid(True)
ax.legend()
plt.show()


# ============================================================================
# 6. MONTE CARLO ERROR VERSUS CORRELATION
# ============================================================================

print_header(6, "MONTE CARLO ERROR VERSUS CORRELATION")

call_max_error = (
    np.array(call_max_mc)
    - np.array(call_max_analytic)
)

call_min_error = (
    np.array(call_min_mc)
    - np.array(call_min_analytic)
)

put_max_error = (
    np.array(put_max_mc)
    - np.array(put_max_analytic)
)

put_min_error = (
    np.array(put_min_mc)
    - np.array(put_min_analytic)
)

fig, ax = plt.subplots()

ax.axhline(
    0.0,
    linestyle="--",
)

ax.plot(
    corr_list,
    call_max_error,
    marker="o",
    label="Call maximum",
)

ax.plot(
    corr_list,
    call_min_error,
    marker="o",
    label="Call minimum",
)

ax.plot(
    corr_list,
    put_max_error,
    marker="o",
    label="Put maximum",
)

ax.plot(
    corr_list,
    put_min_error,
    marker="o",
    label="Put minimum",
)

ax.set_xlabel("Correlation")
ax.set_ylabel("MC - analytic")
ax.set_title("Monte Carlo Error versus Correlation")
ax.grid(True)
ax.legend()
plt.show()


# ============================================================================
# 7. NTH-ASSET CONSISTENCY
# ============================================================================

print_header(7, "NTH-ASSET CONSISTENCY")

correlation = 0.50

corr_matrix = make_corr_matrix(
    num_assets,
    correlation,
)

call_max = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

call_min = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    [strike],
    num_assets,
)

call_first = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_NTH,
    [1, strike],
    num_assets,
)

call_second = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_NTH,
    [2, strike],
    num_assets,
)

value_max = call_max.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
)

value_min = call_min.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
)

value_first_mc = call_first.value_mc(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
    num_paths,
    seed=seed,
)

value_second_mc = call_second.value_mc(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
    num_paths,
    seed=seed,
)

print(f"Correlation       : {correlation:.4f}")
print(f"Paths             : {num_paths}")
print()

print(
    f"{'METHOD':<24}"
    f"{'VALUE':>16}"
    f"{'DIFFERENCE':>16}"
)

print("-" * 56)

print(
    f"{'Call on maximum':<24}"
    f"{value_max:16.8f}"
    f"{0.0:16.8f}"
)

print(
    f"{'Call on 1st MC':<24}"
    f"{value_first_mc:16.8f}"
    f"{value_first_mc - value_max:16.8f}"
)

print(
    f"{'Call on minimum':<24}"
    f"{value_min:16.8f}"
    f"{0.0:16.8f}"
)

print(
    f"{'Call on 2nd MC':<24}"
    f"{value_second_mc:16.8f}"
    f"{value_second_mc - value_min:16.8f}"
)


# ============================================================================
# 8. MONTE CARLO PATH CONVERGENCE
# ============================================================================

print_header(8, "MONTE CARLO PATH CONVERGENCE")

correlation = 0.50

corr_matrix = make_corr_matrix(
    num_assets,
    correlation,
)

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

analytic_value = option.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
)

path_list = [
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
]

num_replications = 20
base_seed = 4242

mc_means = []
mc_stds = []
mc_ses = []
mc_biases = []

print("Option type       : CALL_ON_MAXIMUM")
print(f"Correlation       : {correlation:.4f}")
print(f"Analytic value    : {analytic_value:.10f}")
print(f"Replications      : {num_replications}")
print()

print(
    f"{'PATHS':>10}"
    f"{'MC MEAN':>16}"
    f"{'MC STD':>16}"
    f"{'MEAN SE':>16}"
    f"{'BIAS':>16}"
    f"{'BIAS/SE':>12}"
)

print("-" * 86)

for num_paths_test in path_list:

    values = []

    for replication in range(num_replications):

        replication_seed = base_seed + replication

        value_mc = option.value_mc(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths_test,
            seed=replication_seed,
        )

        values.append(value_mc)

    values = np.array(values)

    mc_mean = np.mean(values)
    mc_std = np.std(values, ddof=1)
    mc_se = mc_std / np.sqrt(num_replications)

    bias = mc_mean - analytic_value

    if mc_se > 0.0:
        bias_se = bias / mc_se
    else:
        bias_se = np.nan

    mc_means.append(mc_mean)
    mc_stds.append(mc_std)
    mc_ses.append(mc_se)
    mc_biases.append(bias)

    print(
        f"{num_paths_test:10d}"
        f"{mc_mean:16.8f}"
        f"{mc_std:16.8f}"
        f"{mc_se:16.8f}"
        f"{bias:16.8f}"
        f"{bias_se:12.4f}"
    )


# ============================================================================
# 9. OPTION VALUE VERSUS NUMBER OF PATHS
# ============================================================================

print_header(9, "OPTION VALUE VERSUS NUMBER OF PATHS")

path_array = np.array(
    path_list,
    dtype=float,
)

mc_mean_array = np.array(mc_means)
mc_se_array = np.array(mc_ses)

fig, ax = plt.subplots()

ax.axhline(
    analytic_value,
    linestyle="--",
    label="Analytic",
)

ax.errorbar(
    path_array,
    mc_mean_array,
    yerr=2.0 * mc_se_array,
    marker="o",
    capsize=4,
    label="MC mean +/- 2 SE",
)

ax.set_xscale("log")

ax.set_xlabel("Number of paths")
ax.set_ylabel("Option value")
ax.set_title("Rainbow Option Monte Carlo Convergence")
ax.grid(True)
ax.legend()
plt.show()


# ============================================================================
# 10. MONTE CARLO ERROR VERSUS NUMBER OF PATHS
# ============================================================================

print_header(
    10,
    "MONTE CARLO ERROR VERSUS NUMBER OF PATHS",
)

mc_bias_array = np.array(mc_biases)

fig, ax = plt.subplots()

ax.axhline(
    0.0,
    linestyle="--",
)

ax.errorbar(
    path_array,
    mc_bias_array,
    yerr=2.0 * mc_se_array,
    marker="o",
    capsize=4,
)

ax.set_xscale("log")

ax.set_xlabel("Number of paths")
ax.set_ylabel("MC mean - analytic")
ax.set_title("Rainbow Option Monte Carlo Error")
ax.grid(True)
plt.show()


# ============================================================================
# 11. MONTE CARLO STATISTICAL CONVERGENCE
# ============================================================================

print_header(
    11,
    "MONTE CARLO STATISTICAL CONVERGENCE",
)

mc_std_array = np.array(mc_stds)

reference_std = (
    mc_std_array[0]
    * np.sqrt(path_array[0] / path_array)
)

fig, ax = plt.subplots()

ax.plot(
    path_array,
    mc_std_array,
    marker="o",
    label="Observed MC std",
)

ax.plot(
    path_array,
    reference_std,
    linestyle="--",
    label=r"$1/\sqrt{N}$ reference",
)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Number of paths")
ax.set_ylabel("Standard deviation")
ax.set_title("Monte Carlo Statistical Convergence")
ax.grid(True)
ax.legend()
plt.show()


# ============================================================================
# 12. FIVE-ASSET NTH-ORDER OPTIONS
# ============================================================================

print_header(
    12,
    "FIVE-ASSET NTH-ORDER OPTIONS",
)

num_assets_5 = 5

stock_prices_5 = (
    np.ones(num_assets_5) * 100.0
)

volatilities_5 = (
    np.ones(num_assets_5) * 0.30
)

dividend_curves_5 = [
    FlatDiscountCurve(
        value_dt,
        dividend_yield,
    )
    for _ in range(num_assets_5)
]

correlation = 0.50

corr_matrix_5 = make_corr_matrix(
    num_assets_5,
    correlation,
)

nth_values = []
nth_call_values = []
nth_put_values = []

print(f"Correlation       : {correlation:.4f}")
print(f"Paths             : {num_paths}")
print()

print(
    f"{'NTH':>8}"
    f"{'CALL MC':>16}"
    f"{'PUT MC':>16}"
)

print("-" * 40)

for n in range(1, num_assets_5 + 1):

    call_option = EquityRainbowOption(
        expiry_dt,
        EquityRainbowOptionTypes.CALL_ON_NTH,
        [n, strike],
        num_assets_5,
    )

    put_option = EquityRainbowOption(
        expiry_dt,
        EquityRainbowOptionTypes.PUT_ON_NTH,
        [n, strike],
        num_assets_5,
    )

    call_value = call_option.value_mc(
        value_dt,
        stock_prices_5,
        discount_curve,
        dividend_curves_5,
        volatilities_5,
        corr_matrix_5,
        num_paths,
        seed=seed,
    )

    put_value = put_option.value_mc(
        value_dt,
        stock_prices_5,
        discount_curve,
        dividend_curves_5,
        volatilities_5,
        corr_matrix_5,
        num_paths,
        seed=seed,
    )

    nth_values.append(n)
    nth_call_values.append(call_value)
    nth_put_values.append(put_value)

    print(
        f"{n:8d}"
        f"{call_value:16.8f}"
        f"{put_value:16.8f}"
    )


# ============================================================================
# 13. OPTION VALUE VERSUS ORDER STATISTIC
# ============================================================================

print_header(
    13,
    "OPTION VALUE VERSUS ORDER STATISTIC",
)

fig, ax = plt.subplots()

ax.plot(
    nth_values,
    nth_call_values,
    marker="o",
    label="Call",
)

ax.plot(
    nth_values,
    nth_put_values,
    marker="o",
    label="Put",
)

ax.set_xlabel("Order statistic")
ax.set_ylabel("Option value")
ax.set_title(
    "Five-Asset Rainbow Option Value versus Order Statistic"
)

ax.set_xticks(nth_values)

ax.grid(True)
ax.legend()
plt.show()


num_paths = 10000
num_replications = 100
base_seed = 4242

option_types = [
    (EquityRainbowOptionTypes.CALL_ON_MAXIMUM, "Call maximum"),
    (EquityRainbowOptionTypes.CALL_ON_MINIMUM, "Call minimum"),
    (EquityRainbowOptionTypes.PUT_ON_MAXIMUM, "Put maximum"),
    (EquityRainbowOptionTypes.PUT_ON_MINIMUM, "Put minimum"),
]

print(
    f"{'TYPE':<18}"
    f"{'RHO':>8}"
    f"{'ANALYTIC':>14}"
    f"{'MC MEAN':>14}"
    f"{'MC STD':>14}"
    f"{'MEAN SE':>14}"
    f"{'BIAS':>14}"
    f"{'BIAS/SE':>12}"
)

print("-" * 108)

for option_type, label in option_types:

    option = EquityRainbowOption(
        expiry_dt,
        option_type,
        [strike],
        num_assets,
    )

    for correlation in corr_list:

        corr_matrix = make_corr_matrix(
            num_assets,
            correlation,
        )

        analytic = option.value(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
        )

        values = []

        for rep in range(num_replications):

            mc = option.value_mc(
                value_dt,
                stock_prices,
                discount_curve,
                dividend_curves,
                volatilities,
                corr_matrix,
                num_paths,
                seed=base_seed + rep,
            )

            values.append(mc)

        values = np.asarray(values)

        mc_mean = np.mean(values)
        mc_std = np.std(values, ddof=1)
        mc_se = mc_std / np.sqrt(num_replications)

        bias = mc_mean - analytic
        bias_se = bias / mc_se

        print(
            f"{label:<18}"
            f"{correlation:8.2f}"
            f"{analytic:14.8f}"
            f"{mc_mean:14.8f}"
            f"{mc_std:14.8f}"
            f"{mc_se:14.8f}"
            f"{bias:14.8f}"
            f"{bias_se:12.4f}"
        )

    print()

# ============================================================================
# ANALYTIC CALL MAX/MIN PARITY
# ============================================================================


print_header(
    14,
    "ANALYTIC CALL MAX/MIN PARITY",
)

call_max = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

call_min = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    [strike],
    num_assets,
)

print(
    f"{'RHO':>8}"
    f"{'CALL MAX':>14}"
    f"{'CALL MIN':>14}"
    f"{'SUM':>14}"
    f"{'VANILLA SUM':>14}"
    f"{'DIFF':>14}"
)

print("-" * 78)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    cmax = call_max.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    cmin = call_min.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    # Black-Scholes vanilla calls on the individual assets.
    vanilla_sum = 0.0

    for i in range(num_assets):

        vanilla_option = EquityVanillaOption(
            expiry_dt,
            strike,
            OptionTypes.EUROPEAN_CALL,
        )

        vanilla_sum += vanilla_option.value(
            value_dt,
            stock_prices[i],
            discount_curve,
            dividend_curves[i],
            BlackScholes(volatilities[i]),
        )

    rainbow_sum = cmax + cmin
    diff = rainbow_sum - vanilla_sum

    print(
        f"{correlation:8.2f}"
        f"{cmax:14.8f}"
        f"{cmin:14.8f}"
        f"{rainbow_sum:14.8f}"
        f"{vanilla_sum:14.8f}"
        f"{diff:14.10f}"
    )


# ============================================================================
# ANALYTIC PUT MAX/MIN PARITY
# ============================================================================

print_header(
    15,
    "ANALYTIC PUT MAX/MIN PARITY",
)

put_max = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MAXIMUM,
    [strike],
    num_assets,
)

put_min = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MINIMUM,
    [strike],
    num_assets,
)

print(
    f"{'RHO':>8}"
    f"{'PUT MAX':>14}"
    f"{'PUT MIN':>14}"
    f"{'SUM':>14}"
    f"{'VANILLA SUM':>14}"
    f"{'DIFF':>14}"
)

print("-" * 78)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    pmax = put_max.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    pmin = put_min.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    vanilla_sum = 0.0

    for i in range(num_assets):

        vanilla_option = EquityVanillaOption(
            expiry_dt,
            strike,
            OptionTypes.EUROPEAN_PUT,
        )

        vanilla_sum += vanilla_option.value(
            value_dt,
            stock_prices[i],
            discount_curve,
            dividend_curves[i],
            BlackScholes(volatilities[i]),
        )

    rainbow_sum = pmax + pmin
    diff = rainbow_sum - vanilla_sum

    print(
        f"{correlation:8.2f}"
        f"{pmax:14.8f}"
        f"{pmin:14.8f}"
        f"{rainbow_sum:14.8f}"
        f"{vanilla_sum:14.8f}"
        f"{diff:14.10f}"
    )


# ============================================================================
# RAINBOW PUT-CALL PARITY
# ============================================================================

print_header(
    16,
    "RAINBOW PUT-CALL PARITY",
)

call_max = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

call_min = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    [strike],
    num_assets,
)

put_max = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MAXIMUM,
    [strike],
    num_assets,
)

put_min = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.PUT_ON_MINIMUM,
    [strike],
    num_assets,
)

t_exp = (expiry_dt - value_dt) / 365.0

r = discount_curve.zero_rate_cc(expiry_dt)

q1 = dividend_curves[0].zero_rate_cc(expiry_dt)
q2 = dividend_curves[1].zero_rate_cc(expiry_dt)

rhs = (
    stock_prices[0] * np.exp(-q1 * t_exp)
    + stock_prices[1] * np.exp(-q2 * t_exp)
    - 2.0 * strike * np.exp(-r * t_exp)
)

print(
    f"{'RHO':>8}"
    f"{'C MAX':>14}"
    f"{'C MIN':>14}"
    f"{'P MAX':>14}"
    f"{'P MIN':>14}"
    f"{'LHS':>14}"
    f"{'RHS':>14}"
    f"{'DIFF':>14}"
)

print("-" * 106)

for correlation in corr_list:

    corr_matrix = make_corr_matrix(
        num_assets,
        correlation,
    )

    cmax = call_max.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    cmin = call_min.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    pmax = put_max.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    pmin = put_min.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    lhs = cmax - pmax + cmin - pmin

    diff = lhs - rhs

    print(
        f"{correlation:8.2f}"
        f"{cmax:14.8f}"
        f"{cmin:14.8f}"
        f"{pmax:14.8f}"
        f"{pmin:14.8f}"
        f"{lhs:14.8f}"
        f"{rhs:14.8f}"
        f"{diff:14.10f}"
    )


# ============================================================================
# 17. CONTROL-VARIATE MONTE CARLO PATH CONVERGENCE
# ============================================================================

print_header(17, "CONTROL-VARIATE MONTE CARLO PATH CONVERGENCE")

correlation = 0.50
corr_matrix = make_corr_matrix(num_assets, correlation)

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

analytic_value = option.value(
    value_dt,
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    corr_matrix,
)

path_list_cv = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
num_replications_cv = 20
base_seed_cv = 4242

mc_means_cv = []
mc_stds_cv = []
mc_ses_cv = []
mc_biases_cv = []
cv_means = []
cv_stds = []
cv_ses = []
cv_biases = []
variance_reductions = []

print("Option type       : CALL_ON_MAXIMUM")
print(f"Correlation       : {correlation:.4f}")
print(f"Analytic value    : {analytic_value:.10f}")
print(f"Replications      : {num_replications_cv}")
print()
print(
    f"{'PATHS':>10}"
    f"{'MC MEAN':>14}"
    f"{'MC STD':>14}"
    f"{'CV MEAN':>14}"
    f"{'CV STD':>14}"
    f"{'MC BIAS':>14}"
    f"{'CV BIAS':>14}"
    f"{'VR':>12}"
)
print("-" * 110)

for num_paths_test in path_list_cv:
    mc_values = []
    cv_values = []

    for replication in range(num_replications_cv):
        replication_seed = base_seed_cv + replication

        mc = option.value_mc(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths_test,
            seed=replication_seed,
        )

        cv = option.value_mc_cv(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths_test,
            seed=replication_seed,
        )

        mc_values.append(mc)
        cv_values.append(cv)

    mc_values = np.asarray(mc_values)
    cv_values = np.asarray(cv_values)

    mc_mean = np.mean(mc_values)
    mc_std = np.std(mc_values, ddof=1)
    mc_se = mc_std / np.sqrt(num_replications_cv)
    mc_bias = mc_mean - analytic_value

    cv_mean = np.mean(cv_values)
    cv_std = np.std(cv_values, ddof=1)
    cv_se = cv_std / np.sqrt(num_replications_cv)
    cv_bias = cv_mean - analytic_value

    vr = (mc_std * mc_std) / (cv_std * cv_std)

    mc_means_cv.append(mc_mean)
    mc_stds_cv.append(mc_std)
    mc_ses_cv.append(mc_se)
    mc_biases_cv.append(mc_bias)
    cv_means.append(cv_mean)
    cv_stds.append(cv_std)
    cv_ses.append(cv_se)
    cv_biases.append(cv_bias)
    variance_reductions.append(vr)

    print(
        f"{num_paths_test:10d}"
        f"{mc_mean:14.8f}"
        f"{mc_std:14.8f}"
        f"{cv_mean:14.8f}"
        f"{cv_std:14.8f}"
        f"{mc_bias:14.8f}"
        f"{cv_bias:14.8f}"
        f"{vr:12.2f}"
    )

path_array_cv = np.asarray(path_list_cv, dtype=float)
mc_mean_array_cv = np.asarray(mc_means_cv)
mc_std_array_cv = np.asarray(mc_stds_cv)
mc_se_array_cv = np.asarray(mc_ses_cv)
mc_bias_array_cv = np.asarray(mc_biases_cv)
cv_mean_array = np.asarray(cv_means)
cv_std_array = np.asarray(cv_stds)
cv_se_array = np.asarray(cv_ses)
cv_bias_array = np.asarray(cv_biases)
vr_array = np.asarray(variance_reductions)

fig, ax = plt.subplots()
ax.axhline(analytic_value, linestyle="--", label="Analytic")
ax.errorbar(
    path_array_cv,
    mc_mean_array_cv,
    yerr=2.0 * mc_se_array_cv,
    marker="o",
    capsize=4,
    label="MC mean +/- 2 SE",
)
ax.errorbar(
    path_array_cv,
    cv_mean_array,
    yerr=2.0 * cv_se_array,
    marker="s",
    capsize=4,
    label="CV MC mean +/- 2 SE",
)
ax.set_xscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Option value")
ax.set_title("Rainbow MC and Control-Variate Convergence")
ax.grid(True)
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.axhline(0.0, linestyle="--")
ax.errorbar(
    path_array_cv,
    mc_bias_array_cv,
    yerr=2.0 * mc_se_array_cv,
    marker="o",
    capsize=4,
    label="MC",
)
ax.errorbar(
    path_array_cv,
    cv_bias_array,
    yerr=2.0 * cv_se_array,
    marker="s",
    capsize=4,
    label="CV MC",
)
ax.set_xscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Mean - analytic")
ax.set_title("Rainbow MC Error versus Number of Paths")
ax.grid(True)
ax.legend()
plt.show()

mc_reference_std = mc_std_array_cv[0] * np.sqrt(path_array_cv[0] / path_array_cv)
cv_reference_std = cv_std_array[0] * np.sqrt(path_array_cv[0] / path_array_cv)

fig, ax = plt.subplots()
ax.plot(path_array_cv, mc_std_array_cv, marker="o", label="MC std")
ax.plot(path_array_cv, cv_std_array, marker="s", label="CV MC std")
ax.plot(path_array_cv, mc_reference_std, linestyle="--", label=r"MC $1/\sqrt{N}$ reference")
ax.plot(path_array_cv, cv_reference_std, linestyle=":", label=r"CV $1/\sqrt{N}$ reference")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Standard deviation")
ax.set_title("MC Statistical Convergence with Control Variates")
ax.grid(True)
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(path_array_cv, vr_array, marker="o")
ax.axhline(1.0, linestyle="--")
ax.set_xscale("log")
ax.set_xlabel("Number of paths")
ax.set_ylabel("Variance reduction factor")
ax.set_title("Control-Variate Variance Reduction")
ax.grid(True)
plt.show()


# ============================================================================
# 18. CONTROL-VARIATE ERROR VERSUS CORRELATION
# ============================================================================

print_header(18, "CONTROL-VARIATE ERROR VERSUS CORRELATION")

num_paths_corr_cv = 10000
num_replications_corr_cv = 50
base_seed_corr_cv = 4242

option = EquityRainbowOption(
    expiry_dt,
    EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
    [strike],
    num_assets,
)

corr_mc_biases = []
corr_cv_biases = []
corr_mc_stds = []
corr_cv_stds = []
corr_vrs = []

print(
    f"{'RHO':>8}"
    f"{'ANALYTIC':>14}"
    f"{'MC MEAN':>14}"
    f"{'CV MEAN':>14}"
    f"{'MC STD':>14}"
    f"{'CV STD':>14}"
    f"{'VR':>12}"
)
print("-" * 90)

for correlation in corr_list:
    corr_matrix = make_corr_matrix(num_assets, correlation)

    analytic = option.value(
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
    )

    mc_values = []
    cv_values = []

    for rep in range(num_replications_corr_cv):
        replication_seed = base_seed_corr_cv + rep

        mc = option.value_mc(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths_corr_cv,
            seed=replication_seed,
        )

        cv = option.value_mc_cv(
            value_dt,
            stock_prices,
            discount_curve,
            dividend_curves,
            volatilities,
            corr_matrix,
            num_paths_corr_cv,
            seed=replication_seed,
        )

        mc_values.append(mc)
        cv_values.append(cv)

    mc_values = np.asarray(mc_values)
    cv_values = np.asarray(cv_values)

    mc_mean = np.mean(mc_values)
    cv_mean = np.mean(cv_values)
    mc_std = np.std(mc_values, ddof=1)
    cv_std = np.std(cv_values, ddof=1)
    vr = (mc_std * mc_std) / (cv_std * cv_std)

    corr_mc_biases.append(mc_mean - analytic)
    corr_cv_biases.append(cv_mean - analytic)
    corr_mc_stds.append(mc_std)
    corr_cv_stds.append(cv_std)
    corr_vrs.append(vr)

    print(
        f"{correlation:8.2f}"
        f"{analytic:14.8f}"
        f"{mc_mean:14.8f}"
        f"{cv_mean:14.8f}"
        f"{mc_std:14.8f}"
        f"{cv_std:14.8f}"
        f"{vr:12.2f}"
    )

fig, ax = plt.subplots()
ax.axhline(0.0, linestyle="--")
ax.plot(corr_list, corr_mc_biases, marker="o", label="MC mean - analytic")
ax.plot(corr_list, corr_cv_biases, marker="s", label="CV MC mean - analytic")
ax.set_xlabel("Correlation")
ax.set_ylabel("Mean error")
ax.set_title("Control-Variate Error versus Correlation")
ax.grid(True)
ax.legend()
plt.show()

fig, ax = plt.subplots()
ax.plot(corr_list, corr_vrs, marker="o")
ax.axhline(1.0, linestyle="--")
ax.set_xlabel("Correlation")
ax.set_ylabel("Variance reduction factor")
ax.set_title("Control-Variate Efficiency versus Correlation")
ax.grid(True)
plt.show()
