# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import time

import numpy as np


from financepy.models.vasicek_mc import zero_price, zero_price_mc
import matplotlib.pyplot as plt

# ============================================================================
# FINANCEPY EXAMPLES - Vasicek Rate Model
# ============================================================================



########################################################################################




########################################################################################

# ============================================================================
# 1. FIN MODEL RATES VASICEK
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. FIN MODEL RATES VASICEK")
print("=" * 78)

r0 = 0.05
a = 0.10
b = 0.05
sigma = 0.05
t = 5.0

p = zero_price(r0, a, b, sigma, t)

num_paths = 1000
dt = 0.02
seed = 1968

print(f"{'Maturity':>10s} {'Analytic':>12s} {'MC 1,000':>12s} {'MC 10,000':>12s} {'Time (s)':>12s}")
print("-" * 64)

for t in np.linspace(0, 10, 21):
    start = time.time()
    p_mc = zero_price_mc(r0, a, b, sigma, t, dt, num_paths, seed)
    p_mc2 = zero_price_mc(r0, a, b, sigma, t, dt, 10 * num_paths, seed)
    p = zero_price(r0, a, b, sigma, t)
    end = time.time()
    elapsed = end - start
    print(f"{t:10.2f} {p:12.8f} {p_mc:12.8f} {p_mc2:12.8f} {elapsed:12.6f}")

# =============================================================================
# 2. COMPARE ANALYTIC AND MONTE CARLO ZERO-COUPON BOND PRICES
# =============================================================================
# The analytic Vasicek price provides a benchmark for Monte Carlo simulation.
# Increasing the number of paths should generally bring the simulated estimate
# closer to the analytic curve, although sampling noise remains.
plot_times = np.linspace(0.0, 10.0, 21)
analytic_prices = [zero_price(r0, a, b, sigma, x) for x in plot_times]
mc_1k = [zero_price_mc(r0, a, b, sigma, x, dt, 1000, seed) for x in plot_times]
mc_10k = [zero_price_mc(r0, a, b, sigma, x, dt, 10000, seed) for x in plot_times]
plt.figure()
plt.plot(plot_times, analytic_prices, label="Analytic")
plt.plot(plot_times, mc_1k, marker="o", label="Monte Carlo: 1,000 paths")
plt.plot(plot_times, mc_10k, marker="x", label="Monte Carlo: 10,000 paths")
plt.xlabel("Maturity (years)")
plt.ylabel("Zero-coupon bond price")
plt.title("Vasicek analytic price versus Monte Carlo")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
