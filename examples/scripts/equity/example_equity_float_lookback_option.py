# ============================================================================
# FINANCEPY EXAMPLES - EquityFloatLookbackOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import numpy as np
import matplotlib.pyplot as plt

from financepy.models.black_scholes import BlackScholes
from financepy.products.equity.equity_float_lookback_option import EquityFloatLookbackOption
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.utils.date import Date
from financepy.utils.global_types import OptionTypes


value_dt = Date(1, 1, 2015)
expiry_dt = Date(1, 1, 2016)
stock_price = 100.0
volatility = 0.30
interest_rate = 0.05
dividend_yield = 0.01

discount_curve = FlatDiscountCurve(value_dt, interest_rate)
dividend_curve = FlatDiscountCurve(value_dt, dividend_yield)
model = BlackScholes(volatility)

steps_per_year_list = [12, 26, 52, 100, 252, 500, 1000, 2000, 4000]
num_paths_list = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
seed = 4242

# ============================================================================
# 1. FLOATING LOOKBACK CALL
# ============================================================================

print("\n" + "=" * 78)
print("1. FLOATING LOOKBACK CALL")
print("=" * 78)

call_option = EquityFloatLookbackOption(expiry_dt, OptionTypes.EUROPEAN_CALL)
stock_min_max = stock_price

analytic_call = call_option.value(
    value_dt, stock_price, discount_curve, dividend_curve,
    model, stock_min_max
)
print(f"Analytic value : {analytic_call:.10f}")

# ============================================================================
# 2. FLOATING CALL - MONTE CARLO PATH CONVERGENCE
# ============================================================================

print("\n" + "=" * 78)
print("2. FLOATING LOOKBACK CALL - MONTE CARLO PATH CONVERGENCE")
print("=" * 78)

num_steps_per_year = 252
call_path_values = []

print(f"{'PATHS':>12s}{'MC VALUE':>16s}{'ANALYTIC':>16s}{'ERROR':>16s}")
print("-" * 60)

for num_paths in num_paths_list:
    value_mc = call_option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve,
        model, stock_min_max, num_paths=num_paths,
        num_steps_per_year=num_steps_per_year, seed=seed
    )
    call_path_values.append(value_mc)
    print(f"{num_paths:12d}{value_mc:16.8f}{analytic_call:16.8f}"
          f"{value_mc - analytic_call:16.8f}")

plt.figure(figsize=(10, 6))
plt.semilogx(num_paths_list, call_path_values, marker="o", label="Monte Carlo")
plt.axhline(analytic_call, linestyle="--", label="Analytic")
plt.xlabel("Number of Monte Carlo Paths")
plt.ylabel("Floating Lookback Call Value")
plt.title("Floating Lookback Call: Monte Carlo Path Convergence")
plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 3. FLOATING CALL - MONITORING-FREQUENCY CONVERGENCE
# ============================================================================

print("\n" + "=" * 78)
print("3. FLOATING LOOKBACK CALL - OBSERVATION-FREQUENCY CONVERGENCE")
print("=" * 78)

num_paths = 20000
call_monitor_values = []

print(f"{'STEPS/YEAR':>12s}{'MC VALUE':>16s}{'ANALYTIC':>16s}{'ERROR':>16s}")
print("-" * 60)

for steps_per_year in steps_per_year_list:
    value_mc = call_option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve,
        model, stock_min_max, num_paths=num_paths,
        num_steps_per_year=steps_per_year, seed=seed
    )
    call_monitor_values.append(value_mc)
    print(f"{steps_per_year:12d}{value_mc:16.8f}{analytic_call:16.8f}"
          f"{value_mc - analytic_call:16.8f}")

plt.figure(figsize=(10, 6))
plt.plot(steps_per_year_list, call_monitor_values, marker="o", label="Monte Carlo")
plt.axhline(analytic_call, linestyle="--", label="Analytic")
plt.xlabel("Observation Steps per Year")
plt.ylabel("Floating Lookback Call Value")
plt.title("Floating Lookback Call: Observation-Frequency Convergence")
plt.grid(True)
plt.legend()
plt.show()

call_errors = np.asarray(call_monitor_values) - analytic_call
plt.figure(figsize=(10, 6))
plt.plot(steps_per_year_list, call_errors, marker="o", label="MC - Analytic")
plt.axhline(0.0, linestyle="--", label="Zero Error")
plt.xlabel("Observation Steps per Year")
plt.ylabel("MC Value - Analytic Value")
plt.title("Floating Lookback Call: Observation Error")
plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 4. FLOATING LOOKBACK PUT
# ============================================================================

print("\n" + "=" * 78)
print("4. FLOATING LOOKBACK PUT")
print("=" * 78)

put_option = EquityFloatLookbackOption(expiry_dt, OptionTypes.EUROPEAN_PUT)
stock_min_max = stock_price

analytic_put = put_option.value(
    value_dt, stock_price, discount_curve, dividend_curve,
    model, stock_min_max
)
print(f"Analytic value : {analytic_put:.10f}")

# ============================================================================
# 5. FLOATING PUT - MONITORING-FREQUENCY CONVERGENCE
# ============================================================================

print("\n" + "=" * 78)
print("5. FLOATING LOOKBACK PUT - OBSERVATION-FREQUENCY CONVERGENCE")
print("=" * 78)

put_monitor_values = []
print(f"{'STEPS/YEAR':>12s}{'MC VALUE':>16s}{'ANALYTIC':>16s}{'ERROR':>16s}")
print("-" * 60)

for steps_per_year in steps_per_year_list:
    value_mc = put_option.value_mc(
        value_dt, stock_price, discount_curve, dividend_curve,
        model, stock_min_max, num_paths=num_paths,
        num_steps_per_year=steps_per_year, seed=seed
    )
    put_monitor_values.append(value_mc)
    print(f"{steps_per_year:12d}{value_mc:16.8f}{analytic_put:16.8f}"
          f"{value_mc - analytic_put:16.8f}")

plt.figure(figsize=(10, 6))
plt.plot(steps_per_year_list, put_monitor_values, marker="o", label="Monte Carlo")
plt.axhline(analytic_put, linestyle="--", label="Analytic")
plt.xlabel("Observations per Year")
plt.ylabel("Floating Lookback Put Value")
plt.title("Floating Lookback Put: Observation-Frequency Convergence")
plt.grid(True)
plt.legend()
plt.show()

put_errors = np.asarray(put_monitor_values) - analytic_put
plt.figure(figsize=(10, 6))
plt.plot(steps_per_year_list, put_errors, marker="o", label="MC - Analytic")
plt.axhline(0.0, linestyle="--", label="Zero Error")
plt.xlabel("Monitoring Steps per Year")
plt.ylabel("MC Value - Analytic Value")
plt.title("Floating Lookback Put: Observation Error")
plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 6. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("6. SUMMARY")
print("=" * 78)
print(f"Floating call analytic : {analytic_call:.10f}")
print(f"Floating call MC       : {call_monitor_values[-1]:.10f}")
print(f"Floating call error    : {call_errors[-1]:.10f}")
print()
print(f"Floating put analytic  : {analytic_put:.10f}")
print(f"Floating put MC        : {put_monitor_values[-1]:.10f}")
print(f"Floating put error     : {put_errors[-1]:.10f}")
