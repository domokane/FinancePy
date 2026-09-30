# ============================================================================
# FINANCEPY EXAMPLES - EquityAsianOption
# ============================================================================
#
# Focused examples comparing the Asian-option approximations and Monte Carlo
# pricers. Monte Carlo methods return an MCResult containing:
#
#     result.value
#     result.std_err
#
# ============================================================================

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.models.black_scholes import BlackScholes
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.products.equity.equity_asian_option import (
    AsianOptionValuationTypes,
    EquityAsianOption,
)
from financepy.utils.global_types import OptionTypes


# ============================================================================
# COMMON MARKET AND CONTRACT DATA
# ============================================================================

value_dt = Date(1, 1, 2014)
start_averaging_dt = Date(1, 6, 2014)
expiry_dt = Date(1, 1, 2015)

stock_price = 100.0
strike_price = 100.0
volatility = 0.20
interest_rate = 0.05
dividend_yield = 0.02

num_obs_per_year = 120
num_paths = 10000
seed = 1991

accrued_average = stock_price * 1.10

model = BlackScholes(volatility)

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

asian_option = EquityAsianOption(
    start_averaging_dt,
    expiry_dt,
    strike_price,
    OptionTypes.EUROPEAN_CALL,
    num_obs_per_year,
)


# ============================================================================
# 1. FIVE-WAY COMPARISON
# ============================================================================
# Kemna-Vorst, Turnbull-Wakeman, Curran, Monte Carlo and Monte Carlo with
# control variate.
# ============================================================================

print("\n" + "=" * 78)
print("1. FIVE-WAY COMPARISON")
print("=" * 78)

value_kv = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.KEMNA_VORST,
    accrued_average,
)

value_tw = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.TURNBULL_WAKEMAN,
    accrued_average,
)

value_curran = asian_option.value(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    AsianOptionValuationTypes.CURRAN,
    accrued_average,
)

result_mc = asian_option.value_mc_fast(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

result_mc_cv = asian_option.value_mc_fast_cv(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

print(
    f"{'METHOD':<25s}"
    f"{'VALUE':>15s}"
    f"{'STD ERROR':>15s}"
)
print("-" * 55)

print(f"{'Kemna-Vorst':<25s}{value_kv:15.8f}{'-':>15s}")
print(f"{'Turnbull-Wakeman':<25s}{value_tw:15.8f}{'-':>15s}")
print(f"{'Curran':<25s}{value_curran:15.8f}{'-':>15s}")
print(
    f"{'Monte Carlo':<25s}"
    f"{result_mc.value:15.8f}"
    f"{result_mc.std_err:15.8f}"
)
print(
    f"{'Monte Carlo CV':<25s}"
    f"{result_mc_cv.value:15.8f}"
    f"{result_mc_cv.std_err:15.8f}"
)


# Direct visual comparison of all five methods.
method_names = [
    "Kemna-Vorst",
    "Turnbull-Wakeman",
    "Curran",
    "Monte Carlo",
    "Monte Carlo CV",
]

method_values = [
    value_kv,
    value_tw,
    value_curran,
    result_mc.value,
    result_mc_cv.value,
]

method_errors = [
    0.0,
    0.0,
    0.0,
    1.96 * result_mc.std_err,
    1.96 * result_mc_cv.std_err,
]

plt.figure(figsize=(10, 6))

plt.errorbar(
    method_names,
    method_values,
    yerr=method_errors,
    marker="o",
    linestyle="none",
    capsize=4,
)

plt.ylabel("Asian Option Value")
plt.title(
    "Kemna-Vorst, Turnbull-Wakeman and Curran "
    "versus Monte Carlo"
)
plt.grid(True, axis="y")
plt.xticks(rotation=20)
plt.show()


# ============================================================================
# 2. CURRAN VERSUS MONTE CARLO CV
# ============================================================================

print("\n" + "=" * 78)
print("2. CURRAN VERSUS MONTE CARLO CV")
print("=" * 78)

difference = value_curran - result_mc_cv.value

print(f"Curran value       : {value_curran:.8f}")
print(f"MC CV value        : {result_mc_cv.value:.8f}")
print(f"MC CV std error    : {result_mc_cv.std_err:.8f}")
print(f"Curran - MC CV     : {difference:.8f}")

if result_mc_cv.std_err > 0.0:
    print(
        f"Difference / MC SE : "
        f"{difference / result_mc_cv.std_err:.4f}"
    )


# ============================================================================
# 3. CURRAN VERSUS MONTE CARLO CV BY NUMBER OF OBSERVATIONS
# ============================================================================
# The averaging period starts on the valuation date. Only the number of
# observations per year changes.
# ============================================================================

print("\n" + "=" * 78)
print("3. CURRAN VERSUS MONTE CARLO CV BY NUMBER OF OBSERVATIONS")
print("=" * 78)

value_dt_obs = Date(1, 1, 2015)
start_averaging_dt_obs = value_dt_obs
expiry_dt_obs = Date(1, 1, 2016)

stock_price_obs = 100.0
strike_price_obs = 100.0

discount_curve_obs = FlatDiscountCurve(
    value_dt_obs,
    0.05,
)

dividend_curve_obs = FlatDiscountCurve(
    value_dt_obs,
    0.01,
)

model_obs = BlackScholes(0.30)

num_obs_list = [
    12,
    26,
    52,
    100,
    252,
    500,
]

curran_obs_values = []
mc_cv_obs_values = []
mc_cv_obs_errors = []

print(
    f"{'OBS/YEAR':>10s}"
    f"{'CURRAN':>15s}"
    f"{'MC CV':>15s}"
    f"{'MC SE':>15s}"
    f"{'DIFFERENCE':>15s}"
)
print("-" * 70)

for observations in num_obs_list:

    option_obs = EquityAsianOption(
        start_averaging_dt_obs,
        expiry_dt_obs,
        strike_price_obs,
        OptionTypes.EUROPEAN_CALL,
        observations,
    )

    curran = option_obs.value(
        value_dt_obs,
        stock_price_obs,
        discount_curve_obs,
        dividend_curve_obs,
        model_obs,
        AsianOptionValuationTypes.CURRAN,
        accrued_average=None,
    )

    result_cv = option_obs.value_mc_fast_cv(
        value_dt_obs,
        stock_price_obs,
        discount_curve_obs,
        dividend_curve_obs,
        model_obs,
        num_paths,
        seed,
        accrued_average=None,
    )

    curran_obs_values.append(curran)
    mc_cv_obs_values.append(result_cv.value)
    mc_cv_obs_errors.append(result_cv.std_err)

    print(
        f"{observations:10d}"
        f"{curran:15.8f}"
        f"{result_cv.value:15.8f}"
        f"{result_cv.std_err:15.8f}"
        f"{curran - result_cv.value:15.8f}"
    )

curran_obs_values = np.asarray(curran_obs_values)
mc_cv_obs_values = np.asarray(mc_cv_obs_values)
mc_cv_obs_errors = np.asarray(mc_cv_obs_errors)

plt.figure(figsize=(9, 6))

plt.plot(
    num_obs_list,
    curran_obs_values,
    marker="o",
    label="Curran",
)

plt.errorbar(
    num_obs_list,
    mc_cv_obs_values,
    yerr=1.96 * mc_cv_obs_errors,
    marker="o",
    capsize=3,
    label="Monte Carlo CV (95% CI)",
)

plt.xscale("log")
plt.xlabel("Observations per Year")
plt.ylabel("Asian Option Value")
plt.title("Curran versus Monte Carlo CV by Observation Frequency")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. CURRAN VERSUS MONTE CARLO CV THROUGH TIME
# ============================================================================
#
# Move the valuation date from the pre-averaging period, through the start of
# averaging, and towards expiry.
#
# The stock price is deliberately held constant at 100 throughout this
# analysis. Once averaging has started, the accrued average is also held at
# 100. This isolates the effect of the passage of time and the shrinking
# remaining averaging period.
#
# ============================================================================

# ============================================================================
# 5. DELTA THROUGH THE LIFE OF THE OPTION
# ============================================================================
#
# Calculate Curran delta from the pre-averaging period, through the start of
# averaging, and towards expiry.
#
# The stock price is held constant at 100. Once averaging has started, the
# accrued average is also held at 100. This isolates the effect of time and
# the shrinking remaining averaging period on delta.
#
# The x-axis is time to expiry in years rather than calendar date.
# ============================================================================

print("\n" + "=" * 78)
print("5. DELTA THROUGH THE LIFE OF THE OPTION")
print("=" * 78)

start_averaging_dt_delta = Date(1, 1, 2015)
expiry_dt_delta = Date(1, 1, 2016)

stock_price_delta = 100.0
strike_price_delta = 100.0
volatility_delta = 0.20
interest_rate_delta = 0.05
dividend_yield_delta = 0.02

option_delta = EquityAsianOption(
    start_averaging_dt_delta,
    expiry_dt_delta,
    strike_price_delta,
    OptionTypes.EUROPEAN_CALL,
    100,
)

model_delta = BlackScholes(volatility_delta)

value_dts_delta = [
    Date(1, 4, 2014),
    Date(1, 7, 2014),
    Date(1, 10, 2014),
    Date(1, 1, 2015),
    Date(1, 3, 2015),
    Date(1, 5, 2015),
    Date(1, 7, 2015),
    Date(1, 9, 2015),
    Date(1, 11, 2015),
    Date(1, 12, 2015),
]

times_to_expiry = []
delta_values = []

print(
    f"{'DATE':>15s}"
    f"{'TIME':>15s}"
    f"{'DELTA':>15s}"
)

print("-" * 45)

for time_value_dt in value_dts_delta:

    time_discount_curve = FlatDiscountCurve(
        time_value_dt,
        interest_rate_delta,
    )

    time_dividend_curve = FlatDiscountCurve(
        time_value_dt,
        dividend_yield_delta,
    )

    if time_value_dt <= start_averaging_dt_delta:
        accrued_average_delta = None
    else:
        accrued_average_delta = stock_price_delta

    delta = option_delta.delta(
        time_value_dt,
        stock_price_delta,
        time_discount_curve,
        time_dividend_curve,
        model_delta,
        method=AsianOptionValuationTypes.CURRAN,
        accrued_average=accrued_average_delta,
    )

    time_to_expiry = (
        expiry_dt_delta - time_value_dt
    ) / 365.0

    times_to_expiry.append(time_to_expiry)
    delta_values.append(delta)

    print(
        f"{str(time_value_dt):>15s}"
        f"{time_to_expiry:15.6f}"
        f"{delta:15.8f}"
    )

times_to_expiry = np.asarray(times_to_expiry)
delta_values = np.asarray(delta_values)

time_to_averaging_start = (
    expiry_dt_delta - start_averaging_dt_delta
) / 365.0

plt.figure(figsize=(10, 6))

plt.plot(
    times_to_expiry,
    delta_values,
    marker="o",
    label="Curran Delta",
)

plt.axvline(
    time_to_averaging_start,
    linestyle="--",
    label="Averaging Starts",
)

plt.xlabel("Time to Expiry (Years)")
plt.ylabel("Delta")
plt.title("Asian Option Delta through Time")

plt.grid(True)
plt.legend()

# Show chronological progression towards expiry:
# large time-to-expiry -> zero
plt.gca().invert_xaxis()

plt.tight_layout()
plt.show()


# ============================================================================
# 5. MONTE CARLO CONVERGENCE
# ============================================================================
# Compare MC and MC CV directly against Kemna-Vorst, Turnbull-Wakeman and
# Curran as the number of Monte Carlo paths increases. The returned standard
# errors provide the MC confidence intervals without repeated simulations.
# ============================================================================

print("\n" + "=" * 78)
print("5. MONTE CARLO CONVERGENCE")
print("=" * 78)

num_paths_list = [
    2000,
    4000,
    8000,
    20000,
    50000,
    100000,
]

mc_values = []
mc_errors = []
mc_cv_values = []
mc_cv_errors = []

print(
    f"{'PATHS':>10s}"
    f"{'MC':>15s}"
    f"{'MC SE':>12s}"
    f"{'MC CV':>15s}"
    f"{'CV SE':>12s}"
)
print("-" * 64)

for paths in num_paths_list:

    result_mc = asian_option.value_mc_fast(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    result_cv = asian_option.value_mc_fast_cv(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        paths,
        seed,
        accrued_average,
    )

    mc_values.append(result_mc.value)
    mc_errors.append(result_mc.std_err)
    mc_cv_values.append(result_cv.value)
    mc_cv_errors.append(result_cv.std_err)

    print(
        f"{paths:10d}"
        f"{result_mc.value:15.8f}"
        f"{result_mc.std_err:12.8f}"
        f"{result_cv.value:15.8f}"
        f"{result_cv.std_err:12.8f}"
    )

mc_values = np.asarray(mc_values)
mc_errors = np.asarray(mc_errors)
mc_cv_values = np.asarray(mc_cv_values)
mc_cv_errors = np.asarray(mc_cv_errors)

plt.figure(figsize=(9, 6))

plt.errorbar(
    num_paths_list,
    mc_values,
    yerr=1.96 * mc_errors,
    marker="o",
    capsize=3,
    label="Monte Carlo (95% CI)",
)

plt.errorbar(
    num_paths_list,
    mc_cv_values,
    yerr=1.96 * mc_cv_errors,
    marker="o",
    capsize=3,
    label="Monte Carlo CV (95% CI)",
)

plt.axhline(
    value_kv,
    linestyle=":",
    label="Kemna-Vorst",
)

plt.axhline(
    value_tw,
    linestyle=":",
    label="Turnbull-Wakeman",
)

plt.axhline(
    value_curran,
    linestyle="--",
    label="Curran",
)

plt.xscale("log")
plt.xlabel("Number of Paths")
plt.ylabel("Asian Option Value")
plt.title("Monte Carlo Convergence")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. CONTROL-VARIATE EFFECTIVENESS
# ============================================================================

print("\n" + "=" * 78)
print("6. CONTROL-VARIATE EFFECTIVENESS")
print("=" * 78)

result_mc = asian_option.value_mc_fast(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

result_mc_cv = asian_option.value_mc_fast_cv(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    model,
    num_paths,
    seed,
    accrued_average,
)

variance_reduction = (
    result_mc.std_err
    / result_mc_cv.std_err
) ** 2

print(f"Paths              : {num_paths}")
print(f"MC value           : {result_mc.value:.10f}")
print(f"MC standard error  : {result_mc.std_err:.10f}")
print(f"MC CV value        : {result_mc_cv.value:.10f}")
print(f"MC CV standard error: {result_mc_cv.std_err:.10f}")
print(f"Variance reduction : {variance_reduction:.2f}x")


# ============================================================================
# 7. STANDARD-ERROR VALIDATION
# ============================================================================
# This is the only example that repeats the calculation across independent
# seeds. It checks the standard error returned by the pricer against the
# empirical dispersion of independent MC estimates.
# ============================================================================

print("\n" + "=" * 78)
print("7. STANDARD-ERROR VALIDATION")
print("=" * 78)

num_replications = 100
seed_start = 1000

values = []
reported_errors = []

for test_seed in range(seed_start, seed_start + num_replications):

    result = asian_option.value_mc_fast_cv(
        value_dt,
        stock_price,
        discount_curve,
        dividend_curve,
        model,
        num_paths,
        test_seed,
        accrued_average,
    )

    values.append(result.value)
    reported_errors.append(result.std_err)

values = np.asarray(values)
reported_errors = np.asarray(reported_errors)

empirical_std = np.std(values, ddof=1)
mean_reported_se = np.mean(reported_errors)

print(f"Paths                  : {num_paths}")
print(f"Replications           : {num_replications}")
print(f"Mean MC value          : {np.mean(values):.10f}")
print(f"Empirical std          : {empirical_std:.10f}")
print(f"Mean reported SE       : {mean_reported_se:.10f}")
print(f"Empirical / reported SE: {empirical_std / mean_reported_se:.4f}")
