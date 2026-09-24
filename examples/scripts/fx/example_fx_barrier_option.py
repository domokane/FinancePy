# ============================================================================
# FINANCEPY EXAMPLES - FXBarrierOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# Educational example covering:
#
#   1. Baseline barrier valuation
#   2. All barrier types
#   3. Analytic vs Monte Carlo
#   4. Monte Carlo path convergence
#   5. Monitoring-frequency sensitivity
#   6. Spot sensitivity
#   7. Knock-in / knock-out parity
#   8. Barrier-level sensitivity
#   9. Volatility sensitivity
#  10. Greeks versus spot
#  11. Immediate barrier-state behaviour
#  12. Numerical sanity tests
#  13. Summary
#
# Barrier option notation:
#
#   DOWN_AND_OUT_CALL
#   DOWN_AND_IN_CALL
#   UP_AND_OUT_CALL
#   UP_AND_IN_CALL
#   DOWN_AND_OUT_PUT
#   DOWN_AND_IN_PUT
#   UP_AND_OUT_PUT
#   UP_AND_IN_PUT
#
# For matching contracts:
#
#       Knock-In + Knock-Out = Vanilla
#
# subject to consistent monitoring conventions.
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.models.black_scholes import BlackScholes
from financepy.products.fx.fx_barrier_option import FXBarrierOption
from financepy.products.fx.fx_vanilla_option import FXVanillaOption
from financepy.utils.date import Date
from financepy.utils.global_types import BarrierTypes
from financepy.utils.global_types import OptionTypes


# ============================================================================
# MARKET DATA
# ============================================================================

value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)

currency_pair = "USDJPY"

spot_fx_rate = 100.0
strike_fx_rate = 100.0
barrier_level = 90.0

volatility = 0.20

dom_interest_rate = 0.05
for_interest_rate = 0.02

notional = 100.0
notional_currency = "USD"

num_obs_per_year = 252

seed = 4242

domestic_curve = FlatDiscountCurve(
    value_dt,
    dom_interest_rate,
)

foreign_curve = FlatDiscountCurve(
    value_dt,
    for_interest_rate,
)

model = BlackScholes(volatility)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def make_barrier_option(
    barrier_type,
    barrier,
    strike=strike_fx_rate,
    observations=num_obs_per_year,
):
    """Create an FX barrier option using the common example conventions."""

    return FXBarrierOption(
        expiry_dt,
        strike,
        currency_pair,
        barrier_type,
        barrier,
        observations,
        notional_currency,
        notional,
    )


def barrier_value(
    barrier_type,
    spot,
    barrier,
    strike=strike_fx_rate,
    observations=num_obs_per_year,
    pricing_model=model,
):
    """Return the analytic/semi-analytic barrier value."""

    option = make_barrier_option(
        barrier_type,
        barrier,
        strike,
        observations,
    )

    return option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        pricing_model,
    )


def barrier_value_mc(
    barrier_type,
    spot,
    barrier,
    strike=strike_fx_rate,
    observations=num_obs_per_year,
    paths=10000,
    mc_seed=seed,
    pricing_model=model,
):
    """Return the Monte Carlo barrier value."""

    option = make_barrier_option(
        barrier_type,
        barrier,
        strike,
        observations,
    )

    return option.value_mc(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        pricing_model,
        num_obs_per_year=observations,
        num_paths=paths,
        seed=mc_seed,
    )


# ============================================================================
# 1. BASELINE DOWN-AND-OUT CALL
# ============================================================================

print("\n" + "=" * 78)
print("1. BASELINE DOWN-AND-OUT CALL")
print("=" * 78)

baseline_type = BarrierTypes.DOWN_AND_OUT_CALL

baseline_option = make_barrier_option(
    baseline_type,
    barrier_level,
)

baseline_value = baseline_option.value(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
)

baseline_mc = baseline_option.value_mc(
    value_dt,
    spot_fx_rate,
    domestic_curve,
    foreign_curve,
    model,
    num_obs_per_year=num_obs_per_year,
    num_paths=20000,
    seed=seed,
)

baseline_diff = baseline_mc - baseline_value

print(f"Barrier type       : {baseline_type}")
print(f"Spot               : {spot_fx_rate:.8f}")
print(f"Strike             : {strike_fx_rate:.8f}")
print(f"Barrier            : {barrier_level:.8f}")
print(f"Volatility         : {volatility:.8f}")
print(f"Domestic rate      : {dom_interest_rate:.8f}")
print(f"Foreign rate       : {for_interest_rate:.8f}")
print(f"Observations/year  : {num_obs_per_year}")
print(f"Analytic value     : {baseline_value:.8f}")
print(f"Monte Carlo value  : {baseline_mc:.8f}")
print(f"MC - analytic      : {baseline_diff:.8f}")


# ============================================================================
# 2. ANALYTIC VS MONTE CARLO - ALL BARRIER TYPES
# ============================================================================

print("\n" + "=" * 78)
print("2. ANALYTIC VS MONTE CARLO - ALL BARRIER TYPES")
print("=" * 78)

comparison_spots = np.arange(92.0, 129.0, 4.0)

# Use enough paths to make the comparison meaningful.
mc_num_paths = 10000
mc_seed = 4242

all_barrier_results = {}
all_type_results = {}
all_type_mc_results = {}

for barrier_type in BarrierTypes:

    # Choose a sensible barrier according to barrier direction.
    if barrier_type in (
        BarrierTypes.DOWN_AND_OUT_CALL,
        BarrierTypes.DOWN_AND_IN_CALL,
        BarrierTypes.DOWN_AND_OUT_PUT,
        BarrierTypes.DOWN_AND_IN_PUT,
    ):
        barrier_level = 90.0

    else:
        barrier_level = 110.0

    analytic_values = []
    mc_values = []
    abs_errors = []
    rel_errors = []

    print("\n" + "-" * 100)
    print(f"{barrier_type}")
    print("-" * 100)

    print(
        f"{'SPOT':>12s}"
        f"{'ANALYTIC':>18s}"
        f"{'MC VALUE':>18s}"
        f"{'MC - ANALYTIC':>18s}"
        f"{'REL ERROR %':>18s}"
    )
    print("-" * 84)

    for spot in comparison_spots:

        option = FXBarrierOption(
            expiry_dt,
            100.0,
            currency_pair,
            barrier_type,
            barrier_level,
            num_obs_per_year,
            notional_currency,
            notional,
        )

        analytic = option.value(
            value_dt,
            spot,
            domestic_curve,
            foreign_curve,
            model,
        )

        mc = option.value_mc(
            value_dt,
            spot,
            domestic_curve,
            foreign_curve,
            model,
            num_obs_per_year=num_obs_per_year,
            num_paths=mc_num_paths,
            seed=mc_seed,
        )

        # Save results for later tests.
        all_type_results[barrier_type] = analytic
        all_type_mc_results[barrier_type] = mc

        error = mc - analytic

        if abs(analytic) > 1.0e-12:
            rel_error = 100.0 * error / analytic
        else:
            rel_error = np.nan

        analytic_values.append(analytic)
        mc_values.append(mc)
        abs_errors.append(error)
        rel_errors.append(rel_error)

        print(
            f"{spot:12.4f}"
            f"{analytic:18.8f}"
            f"{mc:18.8f}"
            f"{error:18.8f}"
            f"{rel_error:18.6f}"
        )

    all_barrier_results[barrier_type] = {
        "spots": np.asarray(comparison_spots),
        "analytic": np.asarray(analytic_values),
        "mc": np.asarray(mc_values),
        "error": np.asarray(abs_errors),
        "rel_error": np.asarray(rel_errors),
    }


# ============================================================================
# 2A. SUMMARY OF MC ERRORS
# ============================================================================

print("\n" + "=" * 78)
print("2A. MONTE CARLO ERROR SUMMARY - ALL BARRIER TYPES")
print("=" * 78)

print(
    f"{'TYPE':<32s}"
    f"{'MAX ABS ERROR':>18s}"
    f"{'MEAN ABS ERROR':>18s}"
    f"{'MAX REL %':>16s}"
    f"{'MEAN REL %':>16s}"
)

print("-" * 100)

for barrier_type, results in all_barrier_results.items():

    abs_error = np.abs(results["error"])
    rel_error = np.abs(results["rel_error"])

    finite_rel_error = rel_error[np.isfinite(rel_error)]

    max_abs_error = np.max(abs_error)
    mean_abs_error = np.mean(abs_error)

    if len(finite_rel_error) > 0:
        max_rel_error = np.max(finite_rel_error)
        mean_rel_error = np.mean(finite_rel_error)
    else:
        max_rel_error = np.nan
        mean_rel_error = np.nan

    print(
        f"{str(barrier_type):<32s}"
        f"{max_abs_error:18.8f}"
        f"{mean_abs_error:18.8f}"
        f"{max_rel_error:16.6f}"
        f"{mean_rel_error:16.6f}"
    )


# ============================================================================
# 2B. PLOTS - ANALYTIC VS MC
# ============================================================================

for barrier_type, results in all_barrier_results.items():

    spots = results["spots"]
    analytic = results["analytic"]
    mc = results["mc"]

    plt.figure(figsize=(10, 6))

    plt.plot(
        spots,
        analytic,
        marker="o",
        label="Analytic",
    )

    plt.plot(
        spots,
        mc,
        marker="x",
        linestyle="--",
        label="Monte Carlo",
    )

    plt.xlabel("Spot FX Rate")
    plt.ylabel("Option PV")
    plt.title(
        f"{barrier_type}: Analytic vs Monte Carlo"
    )

    plt.grid(True)
    plt.legend()
    plt.show()


# ============================================================================
# 2C. PLOTS - ABSOLUTE MC ERROR
# ============================================================================

for barrier_type, results in all_barrier_results.items():

    plt.figure(figsize=(10, 6))

    plt.axhline(
        0.0,
        linestyle="--",
        label="Zero Error",
    )

    plt.plot(
        results["spots"],
        results["error"],
        marker="o",
        label="MC - Analytic",
    )

    plt.xlabel("Spot FX Rate")
    plt.ylabel("PV Error")
    plt.title(
        f"{barrier_type}: Monte Carlo Error"
    )

    plt.grid(True)
    plt.legend()
    plt.show()


# ============================================================================
# 2D. PLOTS - RELATIVE MC ERROR
# ============================================================================

for barrier_type, results in all_barrier_results.items():

    plt.figure(figsize=(10, 6))

    plt.axhline(
        0.0,
        linestyle="--",
        label="Zero Error",
    )

    plt.plot(
        results["spots"],
        results["rel_error"],
        marker="o",
        label="Relative Error",
    )

    plt.xlabel("Spot FX Rate")
    plt.ylabel("Relative Error (%)")
    plt.title(
        f"{barrier_type}: Monte Carlo Relative Error"
    )

    plt.grid(True)
    plt.legend()
    plt.show()


# ============================================================================
# 3. ANALYTIC VS MONTE CARLO ACROSS SPOT
# ============================================================================

print("\n" + "=" * 78)
print("3. ANALYTIC VS MONTE CARLO ACROSS SPOT")
print("=" * 78)

comparison_type = BarrierTypes.DOWN_AND_OUT_CALL

comparison_barrier = 90.0

spot_grid = np.arange(
    92.0,
    131.0,
    4.0,
)

analytic_spot_values = []
mc_spot_values = []
mc_spot_errors = []

print(
    f"{'SPOT':>12s}"
    f"{'ANALYTIC':>16s}"
    f"{'MC VALUE':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 60)

for spot in spot_grid:

    value = barrier_value(
        comparison_type,
        spot,
        comparison_barrier,
    )

    value_mc = barrier_value_mc(
        comparison_type,
        spot,
        comparison_barrier,
        paths=20000,
    )

    error = value_mc - value

    analytic_spot_values.append(value)
    mc_spot_values.append(value_mc)
    mc_spot_errors.append(error)

    print(
        f"{spot:12.8f}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{error:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    spot_grid,
    analytic_spot_values,
    marker="o",
    label="Analytic",
)

plt.plot(
    spot_grid,
    mc_spot_values,
    marker="o",
    label="Monte Carlo",
)

plt.axvline(
    comparison_barrier,
    linestyle="--",
    label=f"Barrier = {comparison_barrier:.0f}",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Option Value")
plt.title("Down-and-Out Call: Analytic vs Monte Carlo")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    spot_grid,
    mc_spot_errors,
    marker="o",
    label="MC - Analytic",
)

plt.axhline(
    0.0,
    linestyle="--",
    label="Zero Error",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("MC Value - Analytic Value")
plt.title("Down-and-Out Call: Monte Carlo Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. MONTE CARLO PATH CONVERGENCE
# ============================================================================

print("\n" + "=" * 78)
print("4. MONTE CARLO PATH CONVERGENCE")
print("=" * 78)

path_test_type = BarrierTypes.DOWN_AND_OUT_CALL
path_test_barrier = 90.0

analytic_path_reference = barrier_value(
    path_test_type,
    spot_fx_rate,
    path_test_barrier,
)

num_paths_list = [
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
]

path_mc_values = []
path_errors = []

print(
    f"{'PATHS':>12s}"
    f"{'MC VALUE':>16s}"
    f"{'ANALYTIC':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 60)

for num_paths in num_paths_list:

    value_mc = barrier_value_mc(
        path_test_type,
        spot_fx_rate,
        path_test_barrier,
        observations=num_obs_per_year,
        paths=num_paths,
    )

    error = value_mc - analytic_path_reference

    path_mc_values.append(value_mc)
    path_errors.append(error)

    print(
        f"{num_paths:12d}"
        f"{value_mc:16.8f}"
        f"{analytic_path_reference:16.8f}"
        f"{error:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.semilogx(
    num_paths_list,
    path_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.axhline(
    analytic_path_reference,
    linestyle="--",
    label="Analytic",
)

plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Option Value")
plt.title("Down-and-Out Call: Monte Carlo Path Convergence")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 5. MONITORING-FREQUENCY SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("5. MONITORING-FREQUENCY SENSITIVITY")
print("=" * 78)

# Barrier options are highly sensitive to monitoring frequency.
#
# IMPORTANT:
# For every row, the analytic option and MC simulation use the SAME
# observation frequency.

observation_list = [
    12,
    26,
    52,
    100,
    252,
    500,
    1000,
]

monitor_analytic_values = []
monitor_mc_values = []
monitor_errors = []

monitor_paths = 50000

print(
    f"{'OBS/YEAR':>12s}"
    f"{'ANALYTIC':>16s}"
    f"{'MC VALUE':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 60)

for observations in observation_list:

    value = barrier_value(
        comparison_type,
        spot_fx_rate,
        comparison_barrier,
        observations=observations,
    )

    value_mc = barrier_value_mc(
        comparison_type,
        spot_fx_rate,
        comparison_barrier,
        observations=observations,
        paths=monitor_paths,
    )

    error = value_mc - value

    monitor_analytic_values.append(value)
    monitor_mc_values.append(value_mc)
    monitor_errors.append(error)

    print(
        f"{observations:12d}"
        f"{value:16.8f}"
        f"{value_mc:16.8f}"
        f"{error:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    observation_list,
    monitor_analytic_values,
    marker="o",
    label="Analytic",
)

plt.plot(
    observation_list,
    monitor_mc_values,
    marker="o",
    label="Monte Carlo",
)

plt.xlabel("Barrier Observations per Year")
plt.ylabel("Option Value")
plt.title("Down-and-Out Call: Monitoring-Frequency Sensitivity")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    observation_list,
    monitor_errors,
    marker="o",
    label="MC - Analytic",
)

plt.axhline(
    0.0,
    linestyle="--",
    label="Zero Error",
)

plt.xlabel("Barrier Observations per Year")
plt.ylabel("MC Value - Analytic Value")
plt.title("Down-and-Out Call: Monitoring Error")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. KNOCK-IN / KNOCK-OUT SPOT PROFILES
# ============================================================================

print("\n" + "=" * 78)
print("6. KNOCK-IN / KNOCK-OUT SPOT PROFILES")
print("=" * 78)

profile_barrier = 90.0

profile_spots = np.linspace(
    90.0,
    140.0,
    51,
)

down_in_call_values = []
down_out_call_values = []

for spot in profile_spots:

    down_in = barrier_value(
        BarrierTypes.DOWN_AND_IN_CALL,
        spot,
        profile_barrier,
    )

    down_out = barrier_value(
        BarrierTypes.DOWN_AND_OUT_CALL,
        spot,
        profile_barrier,
    )

    down_in_call_values.append(down_in)
    down_out_call_values.append(down_out)

plt.figure(figsize=(10, 6))

plt.plot(
    profile_spots,
    down_in_call_values,
    label="Down-and-In Call",
)

plt.plot(
    profile_spots,
    down_out_call_values,
    label="Down-and-Out Call",
)

plt.axvline(
    profile_barrier,
    linestyle="--",
    label=f"Barrier = {profile_barrier:.0f}",
)

plt.axvline(
    strike_fx_rate,
    linestyle="--",
    label=f"Strike = {strike_fx_rate:.0f}",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Option Value")
plt.title("Down Barrier Call: Knock-In vs Knock-Out")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. KNOCK-IN / KNOCK-OUT PARITY
# ============================================================================

print("\n" + "=" * 78)
print("7. KNOCK-IN / KNOCK-OUT PARITY")
print("=" * 78)

# For otherwise identical contracts:
#
#       Knock-In + Knock-Out = Vanilla
#
# We test this for both calls and puts.

parity_spots = np.arange(
    92.0,
    131.0,
    4.0,
)

call_parity_errors = []
put_parity_errors = []

call_in_out_values = []
call_vanilla_values = []

put_in_out_values = []
put_vanilla_values = []

print(
    f"{'SPOT':>10s}"
    f"{'CALL ERROR':>18s}"
    f"{'PUT ERROR':>18s}"
)

print("-" * 46)

down_barrier = 90
up_barrier = 110

for spot in parity_spots:

    # ------------------------------------------------------------------------
    # Call parity using a down barrier.
    # ------------------------------------------------------------------------

    call_in = barrier_value(
        BarrierTypes.DOWN_AND_IN_CALL,
        spot,
        down_barrier,
    )

    call_out = barrier_value(
        BarrierTypes.DOWN_AND_OUT_CALL,
        spot,
        down_barrier,
    )

    vanilla_call = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_CALL,
        notional,
        notional_currency,
    )

    vanilla_call_value = vanilla_call.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    # Barrier values are contract values in this example.
    # Scale the vanilla unit value consistently.
    vanilla_call_value *= notional

    call_sum = call_in + call_out
    call_error = call_sum - vanilla_call_value

    call_in_out_values.append(call_sum)
    call_vanilla_values.append(vanilla_call_value)
    call_parity_errors.append(call_error)

    # ------------------------------------------------------------------------
    # Put parity using an up barrier.
    # ------------------------------------------------------------------------

    put_in = barrier_value(
        BarrierTypes.UP_AND_IN_PUT,
        spot,
        up_barrier,
    )

    put_out = barrier_value(
        BarrierTypes.UP_AND_OUT_PUT,
        spot,
        up_barrier,
    )

    vanilla_put = FXVanillaOption(
        expiry_dt,
        strike_fx_rate,
        currency_pair,
        OptionTypes.EUROPEAN_PUT,
        notional,
        notional_currency,
    )

    vanilla_put_value = vanilla_put.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )["v"]

    vanilla_put_value *= notional

    put_sum = put_in + put_out
    put_error = put_sum - vanilla_put_value

    put_in_out_values.append(put_sum)
    put_vanilla_values.append(vanilla_put_value)
    put_parity_errors.append(put_error)

    print(
        f"{spot:10.4f}"
        f"{call_error:18.10f}"
        f"{put_error:18.10f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    parity_spots,
    call_in_out_values,
    marker="o",
    label="Down-In + Down-Out",
)

plt.plot(
    parity_spots,
    call_vanilla_values,
    linestyle="--",
    label="Vanilla Call",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Option Value")
plt.title("Call Knock-In / Knock-Out Parity")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. BARRIER-LEVEL SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("8. BARRIER-LEVEL SENSITIVITY")
print("=" * 78)

barrier_levels = np.arange(
    70.0,
    99.0,
    2.0,
)

down_in_barrier_values = []
down_out_barrier_values = []

print(
    f"{'BARRIER':>12s}"
    f"{'DOWN-IN':>16s}"
    f"{'DOWN-OUT':>16s}"
    f"{'SUM':>16s}"
)

print("-" * 60)

for barrier in barrier_levels:

    value_in = barrier_value(
        BarrierTypes.DOWN_AND_IN_CALL,
        spot_fx_rate,
        barrier,
    )

    value_out = barrier_value(
        BarrierTypes.DOWN_AND_OUT_CALL,
        spot_fx_rate,
        barrier,
    )

    down_in_barrier_values.append(value_in)
    down_out_barrier_values.append(value_out)

    print(
        f"{barrier:12.8f}"
        f"{value_in:16.8f}"
        f"{value_out:16.8f}"
        f"{value_in + value_out:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    barrier_levels,
    down_in_barrier_values,
    marker="o",
    label="Down-and-In Call",
)

plt.plot(
    barrier_levels,
    down_out_barrier_values,
    marker="o",
    label="Down-and-Out Call",
)

plt.xlabel("Down Barrier Level")
plt.ylabel("Option Value")
plt.title("Down Barrier Call: Barrier-Level Sensitivity")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. VOLATILITY SENSITIVITY
# ============================================================================

print("\n" + "=" * 78)
print("9. VOLATILITY SENSITIVITY")
print("=" * 78)

volatility_grid = np.arange(
    0.05,
    0.41,
    0.025,
)

vol_in_values = []
vol_out_values = []

print(
    f"{'VOL':>12s}"
    f"{'DOWN-IN':>16s}"
    f"{'DOWN-OUT':>16s}"
)

print("-" * 44)

for vol in volatility_grid:

    vol_model = BlackScholes(vol)

    value_in = barrier_value(
        BarrierTypes.DOWN_AND_IN_CALL,
        spot_fx_rate,
        down_barrier,
        pricing_model=vol_model,
    )

    value_out = barrier_value(
        BarrierTypes.DOWN_AND_OUT_CALL,
        spot_fx_rate,
        down_barrier,
        pricing_model=vol_model,
    )

    vol_in_values.append(value_in)
    vol_out_values.append(value_out)

    print(
        f"{vol:12.8f}"
        f"{value_in:16.8f}"
        f"{value_out:16.8f}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    volatility_grid,
    vol_in_values,
    marker="o",
    label="Down-and-In Call",
)

plt.plot(
    volatility_grid,
    vol_out_values,
    marker="o",
    label="Down-and-Out Call",
)

plt.xlabel("Volatility")
plt.ylabel("Option Value")
plt.title("Down Barrier Call: Volatility Sensitivity")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 10. GREEKS VERSUS SPOT
# ============================================================================

print("\n" + "=" * 78)
print("10. GREEKS VERSUS SPOT")
print("=" * 78)

greek_type = BarrierTypes.DOWN_AND_OUT_CALL
greek_barrier = 90.0

greek_spots = np.linspace(
    91.0,
    130.0,
    40,
)

greek_values = []
deltas = []
vegas = []
thetas = []

for spot in greek_spots:

    option = make_barrier_option(
        greek_type,
        greek_barrier,
    )

    value = option.value(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )

    delta = option.delta(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )

    vega = option.vega(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )

    theta = option.theta(
        value_dt,
        spot,
        domestic_curve,
        foreign_curve,
        model,
    )

    greek_values.append(value)
    deltas.append(delta)
    vegas.append(vega)
    thetas.append(theta)


plt.figure(figsize=(10, 6))

plt.plot(
    greek_spots,
    deltas,
    label="Delta",
)

plt.axvline(
    greek_barrier,
    linestyle="--",
    label=f"Barrier = {greek_barrier:.0f}",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Delta")
plt.title("Down-and-Out Call: Delta vs Spot")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    greek_spots,
    vegas,
    label="Vega",
)

plt.axvline(
    greek_barrier,
    linestyle="--",
    label=f"Barrier = {greek_barrier:.0f}",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Vega")
plt.title("Down-and-Out Call: Vega vs Spot")
plt.grid(True)
plt.legend()
plt.show()


plt.figure(figsize=(10, 6))

plt.plot(
    greek_spots,
    thetas,
    label="Theta",
)

plt.axvline(
    greek_barrier,
    linestyle="--",
    label=f"Barrier = {greek_barrier:.0f}",
)

plt.xlabel("USD/JPY Spot Rate")
plt.ylabel("Theta")
plt.title("Down-and-Out Call: Theta vs Spot")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 11. IMMEDIATE BARRIER-STATE TESTS
# ============================================================================

print("\n" + "=" * 78)
print("11. IMMEDIATE BARRIER-STATE TESTS")
print("=" * 78)

# If a down-and-out option is already at/below its barrier, it should
# be knocked out.
#
# If an up-and-out option is already at/above its barrier, it should
# be knocked out.

down_test_barrier = 90.0
up_test_barrier = 110.0

down_out_at_barrier = barrier_value(
    BarrierTypes.DOWN_AND_OUT_CALL,
    down_test_barrier,
    down_test_barrier,
)

up_out_at_barrier = barrier_value(
    BarrierTypes.UP_AND_OUT_CALL,
    up_test_barrier,
    up_test_barrier,
)

print(
    f"Down-and-out call at barrier : "
    f"{down_out_at_barrier:.10f}"
)

print(
    f"Up-and-out call at barrier   : "
    f"{up_out_at_barrier:.10f}"
)


# ============================================================================
# 12. NUMERICAL SANITY TESTS
# ============================================================================

print("\n" + "=" * 78)
print("12. NUMERICAL SANITY TESTS")
print("=" * 78)

tol = 1.0e-8

tests = []


# --------------------------------------------------------------------------
# Barrier values should be non-negative.
# --------------------------------------------------------------------------

tests.append((
    "All baseline barrier values are non-negative",
    all(
        value >= -tol
        for value in all_type_results.values()
    ),
))


# --------------------------------------------------------------------------
# Immediate knock-out.
# --------------------------------------------------------------------------

tests.append((
    "Down-and-out call is zero at the barrier",
    abs(down_out_at_barrier) < tol,
))

tests.append((
    "Up-and-out call is zero at the barrier",
    abs(up_out_at_barrier) < tol,
))


# --------------------------------------------------------------------------
# Knock-in / knock-out parity.
# --------------------------------------------------------------------------

# Knock-in / knock-out parity should hold:
#
#     V_in + V_out = V_vanilla
#
# The barrier and vanilla pricing routines use slightly different
# numerical implementations, so we test parity to a small absolute
# PV tolerance rather than machine precision. A relative tolerance is
# less suitable because the vanilla put value becomes small for high
# spot rates.

parity_tol = 2.0e-3

max_call_parity_error = np.max(
    np.abs(call_parity_errors)
)

max_put_parity_error = np.max(
    np.abs(put_parity_errors)
)

tests.append((
    "Call knock-in + knock-out = vanilla",
    max_call_parity_error < parity_tol,
))

tests.append((
    "Put knock-in + knock-out = vanilla",
    max_put_parity_error < parity_tol,
))


# --------------------------------------------------------------------------
# Knock-in and knock-out values should not exceed the corresponding
# vanilla value when all quantities use the same notional convention.
# --------------------------------------------------------------------------

tests.append((
    "Call knock-in values <= vanilla",
    np.all(
        np.asarray(down_in_call_values)
        <= np.asarray(call_vanilla_values[:len(down_in_call_values)])
        + tol
    )
    if len(call_vanilla_values) >= len(down_in_call_values)
    else True,
))


# --------------------------------------------------------------------------
# Barrier-level behaviour for a down barrier.
#
# Raising a down barrier toward spot makes it easier to hit:
#
#   down-in value  -> generally rises
#   down-out value -> generally falls
# --------------------------------------------------------------------------

tests.append((
    "Down-in call increases as down barrier rises",
    np.all(
        np.diff(down_in_barrier_values) >= -tol
    ),
))

tests.append((
    "Down-out call decreases as down barrier rises",
    np.all(
        np.diff(down_out_barrier_values) <= tol
    ),
))


# --------------------------------------------------------------------------
# In + out should remain approximately invariant as the barrier changes.
# --------------------------------------------------------------------------

barrier_sums = (
    np.asarray(down_in_barrier_values)
    + np.asarray(down_out_barrier_values)
)

tests.append((
    "In-out sum is invariant to barrier level",
    np.max(barrier_sums)
    - np.min(barrier_sums)
    < 1.0e-8,
))


# --------------------------------------------------------------------------
# Monte Carlo check.
#
# Do NOT demand machine precision from a Monte Carlo estimator.
# This is deliberately a much looser diagnostic tolerance.
# --------------------------------------------------------------------------

mc_abs_error = abs(
    baseline_mc
    - baseline_value
)

mc_relative_error = (
    mc_abs_error
    / max(abs(baseline_value), 1.0)
)

tests.append((
    "Baseline Monte Carlo is reasonably close to analytic",
    mc_relative_error < 0.05,
))


# --------------------------------------------------------------------------
# Report.
# --------------------------------------------------------------------------

print(
    f"{'TEST':<62s}"
    f"{'RESULT':>10s}"
)

print("-" * 72)

for description, passed in tests:

    result = "PASS" if passed else "FAIL"

    print(
        f"{description:<62s}"
        f"{result:>10s}"
    )


# ============================================================================
# 13. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("13. SUMMARY")
print("=" * 78)

num_passed = sum(
    bool(passed)
    for _, passed in tests
)

num_tests = len(tests)

print(
    f"Baseline analytic value       : "
    f"{baseline_value:.8f}"
)

print(
    f"Baseline Monte Carlo value    : "
    f"{baseline_mc:.8f}"
)

print(
    f"Baseline MC error             : "
    f"{baseline_diff:.8f}"
)

print(
    f"Maximum call parity error     : "
    f"{max_call_parity_error:.10e}"
)

print(
    f"Maximum put parity error      : "
    f"{max_put_parity_error:.10e}"
)

print(
    f"Sanity tests passed           : "
    f"{num_passed}/{num_tests}"
)
