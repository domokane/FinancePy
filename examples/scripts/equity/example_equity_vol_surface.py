# ============================================================================
# FINANCEPY EXAMPLES - EquityVolSurface
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the construction and analysis of an equity
# implied-volatility surface using FinancePy.
#
# The market data consist of implied volatilities observed across:
#
#     - option strikes
#     - option expiries
#
# An SVI volatility function is fitted to each expiry slice.
#
# The example demonstrates:
#
#   1. Market volatility data
#   2. Input volatility smiles
#   3. Construction of an SVI volatility surface
#   4. Calibrated volatility curves
#   5. Volatility and strike as functions of delta
#   6. Volatility curves in delta space
#   7. Implied risk-neutral densities
#   8. Risk-neutral density plots
#   9. Density normalisation checks
#  10. Summary and timing information
#
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils import Date
from financepy.utils.global_types import VolFuncTypes

from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.volatility.equity_vol_surface import EquityVolSurface

# ============================================================================
# 1. MARKET DATA
# ============================================================================
#
# The volatility surface below contains market implied volatilities for
# several strikes and expiries.
#
# Rows correspond to expiries.
# Columns correspond to strikes.
#
# The original volatility data are expressed in percent and are converted
# below into decimal volatility units.
#
# ============================================================================

print("\n" + "=" * 78)
print("1. MARKET DATA")
print("=" * 78)

value_dt = Date(11, 1, 2021)

stock_price = 3800.0

expiry_dts = [
    Date(11, 2, 2021),
    Date(11, 3, 2021),
    Date(11, 4, 2021),
    Date(11, 7, 2021),
    Date(11, 10, 2021),
    Date(11, 1, 2022),
    Date(11, 1, 2023),
]

strikes = np.array(
    [
        3037.0,
        3418.0,
        3608.0,
        3703.0,
        3798.0,
        3893.0,
        3988.0,
        4178.0,
        4557.0,
    ]
)

vol_surface_pct = np.array(
    [
        [42.94, 31.30, 25.88, 22.94, 19.72, 16.90, 15.31, 17.54, 25.67],
        [37.01, 28.25, 24.19, 21.93, 19.57, 17.45, 15.89, 15.34, 21.15],
        [34.68, 27.38, 23.82, 21.85, 19.83, 17.98, 16.52, 15.31, 18.94],
        [31.41, 26.25, 23.51, 22.05, 20.61, 19.25, 18.03, 16.01, 15.90],
        [29.91, 25.58, 23.21, 22.01, 20.83, 19.70, 18.62, 16.63, 14.94],
        [29.26, 25.24, 23.03, 21.91, 20.81, 19.73, 18.69, 16.76, 14.63],
        [27.59, 24.33, 22.72, 21.93, 21.17, 20.43, 19.71, 18.36, 16.26],
    ]
)

# FinancePy expects volatility in decimal form rather than percent.
vol_surface = vol_surface_pct / 100.0


# ============================================================================
# INTEREST-RATE AND DIVIDEND CURVES
# ============================================================================
#
# For simplicity, this example assumes constant interest and dividend rates.
#
# The discount curve is used to discount option cashflows.
#
# The dividend curve represents the continuous dividend yield of the equity
# index.
#
# ============================================================================

interest_rate = 0.020
dividend_yield = 0.010

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

print(f"Value date       : {value_dt}")
print(f"Stock price      : {stock_price:.2f}")
print(f"Interest rate    : {interest_rate:.4%}")
print(f"Dividend yield   : {dividend_yield:.4%}")
print(f"Number expiries  : {len(expiry_dts)}")
print(f"Number strikes   : {len(strikes)}")


# ============================================================================
# 2. INPUT MARKET VOLATILITY SMILES
# ============================================================================
#
# Before calibrating a volatility model it is useful to inspect the raw
# market volatility data.
#
# Each line represents the implied-volatility smile for one expiry.
#
# The vertical line marks the current stock price.
#
# Short-dated equity-index options often exhibit substantial skew, with
# lower-strike options trading at considerably higher implied volatility.
#
# ============================================================================

print("\n" + "=" * 78)
print("2. INPUT MARKET VOLATILITY SMILES")
print("=" * 78)

plt.figure(figsize=(10, 6))

for i, expiry_dt in enumerate(expiry_dts):

    plt.plot(
        strikes,
        100.0 * vol_surface[i],
        marker="o",
        label=str(expiry_dt),
    )

plt.axvline(
    stock_price,
    linestyle="--",
    label="Spot",
)

plt.xlabel("Strike")
plt.ylabel("Implied Volatility (%)")
plt.title("Equity Market Volatility Smiles")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. BUILD THE SVI VOLATILITY SURFACE
# ============================================================================
#
# FinancePy supports several functional forms for fitting volatility smiles.
#
# In this example we use SVI:
#
#       Stochastic Volatility Inspired
#
# SVI is a commonly used parameterisation of equity volatility smiles because
# it can represent both skew and smile curvature.
#
# EquityVolSurface calibrates the chosen volatility function to each expiry
# slice of the market volatility data.
#
# ============================================================================

print("\n" + "=" * 78)
print("3. BUILD SVI VOLATILITY SURFACE")
print("=" * 78)

vol_function_type = VolFuncTypes.SVI

start_time = time.time()

equity_surface = EquityVolSurface(
    value_dt,
    stock_price,
    discount_curve,
    dividend_curve,
    expiry_dts,
    strikes,
    vol_surface,
    vol_function_type,
)

calibration_time = time.time() - start_time

print(f"Volatility function : {vol_function_type}")
print(f"Calibration time    : {calibration_time:.6f} seconds")


# ============================================================================
# 4. CALIBRATED VOLATILITY CURVES
# ============================================================================
#
# EquityVolSurface provides a convenience function for plotting the fitted
# volatility curves.
#
# The plot allows the fitted SVI curves to be compared visually with the
# original market volatility observations.
#
# ============================================================================

print("\n" + "=" * 78)
print("4. CALIBRATED VOLATILITY CURVES")
print("=" * 78)

equity_surface.plot_vol_curves()


# ============================================================================
# 5. VOLATILITY FROM DELTA
# ============================================================================
#
# Equity-option volatility is often quoted in terms of option delta rather
# than strike.
#
# Given:
#
#       delta
#       expiry
#
# the volatility surface can be inverted to determine both:
#
#       implied volatility
#       corresponding strike
#
# We calculate these quantities for deltas from 10% to 90%.
#
# ============================================================================

print("\n" + "=" * 78)
print("5. VOLATILITY FROM DELTA")
print("=" * 78)

deltas = np.linspace(
    0.10,
    0.90,
    9,
)

print(
    f"{'EXPIRY':>14s}"
    f"{'DELTA':>10s}"
    f"{'VOL':>14s}"
    f"{'STRIKE':>14s}"
)

print("-" * 52)

for expiry_dt in expiry_dts:

    for delta in deltas:

        vol, strike = equity_surface.vol_from_delta_date(
            delta,
            expiry_dt,
        )

        print(
            f"{str(expiry_dt):>14s}"
            f"{delta:10.2f}"
            f"{vol:14.6f}"
            f"{strike:14.4f}"
        )


# ============================================================================
# 6. VOLATILITY CURVES IN DELTA SPACE
# ============================================================================
#
# The volatility surface can also be interrogated using option delta rather
# than strike.
#
# For each expiry and delta, vol_from_delta_date() determines the strike
# corresponding to that delta and evaluates the fitted volatility smile at
# that strike.
#
# Extreme deltas can correspond to strikes far outside the range of the
# market observations. The resulting volatility then depends heavily on
# extrapolation of the fitted SVI smile.
#
# We therefore use deltas between 20% and 80%, which keeps the example
# primarily within the interpolation region while still allowing modest
# extrapolation at the wings.
#
# ============================================================================

print("\n" + "=" * 78)
print("6. VOLATILITY CURVES IN DELTA SPACE")
print("=" * 78)

deltas = np.linspace(0.20, 0.80, 7)

print(
    f"{'EXPIRY':>14s}"
    f"{'DELTA':>10s}"
    f"{'VOL (%)':>14s}"
    f"{'STRIKE':>14s}"
)

print("-" * 52)

plt.figure(figsize=(10, 6))

for expiry_dt in expiry_dts:

    delta_vols = []

    for delta in deltas:

        vol, strike = equity_surface.vol_from_delta_date(
            delta,
            expiry_dt,
        )

        delta_vols.append(vol)

        print(
            f"{str(expiry_dt):>14s}"
            f"{delta:10.2f}"
            f"{100.0 * vol:14.4f}"
            f"{strike:14.2f}"
        )

    plt.plot(
        deltas,
        100.0 * np.asarray(delta_vols),
        marker="o",
        label=str(expiry_dt),
    )

plt.xlabel("Delta")
plt.ylabel("Implied Volatility (%)")
plt.title("Equity Volatility Surface in Delta Space")
plt.grid(True)
plt.legend()
plt.show()

# ============================================================================
# 7. IMPLIED RISK-NEUTRAL DENSITIES
# ============================================================================
#
# European option prices across strike contain information about the
# risk-neutral distribution of the underlying asset at expiry.
#
# FinancePy can calculate the implied probability density associated with
# each expiry slice of the calibrated volatility surface.
#
# The density is evaluated over a strike range substantially wider than the
# original market data.
#
# A sufficiently wide range is important because the probability density
# should integrate to approximately one.
#
# ============================================================================

print("\n" + "=" * 78)
print("7. IMPLIED RISK-NEUTRAL DENSITIES")
print("=" * 78)

min_strike = strikes[0] * 0.5
max_strike = strikes[-1] * 1.5

num_density_points = 1000

dbns = equity_surface.implied_dbns(
    min_strike,
    max_strike,
    num_density_points,
)

print(
    f"{'EXPIRY':>14s}"
    f"{'DENSITY SUM':>18s}"
)

print("-" * 32)

for i, dbn in enumerate(dbns):

    print(
        f"{str(expiry_dts[i]):>14s}"
        f"{dbn.sum():18.10f}"
    )


# ============================================================================
# 8. PLOT IMPLIED RISK-NEUTRAL DENSITIES
# ============================================================================
#
# The risk-neutral densities provide a different view of the information
# contained in the volatility surface.
#
# Differences in skew and curvature translate into differences in the shape
# of the implied probability distributions.
#
# ============================================================================

print("\n" + "=" * 78)
print("8. IMPLIED RISK-NEUTRAL DENSITY PLOTS")
print("=" * 78)

plt.figure(figsize=(10, 6))

for i, dbn in enumerate(dbns):

    plt.plot(
        dbn._x,
        dbn._densitydx,
        label=str(expiry_dts[i]),
    )

plt.xlabel("Stock Price at Expiry")
plt.ylabel("Probability Density")
plt.title("Risk-Neutral Densities Implied by SVI Surface")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 9. DENSITY NORMALISATION CHECK
# ============================================================================
#
# A probability density should integrate to one:
#
#                     integral p(S) dS = 1
#
# The Distribution.sum() method provides a numerical approximation to this
# integral.
#
# The difference
#
#                     density sum - 1
#
# therefore provides a useful numerical diagnostic.
#
# Small errors can arise because:
#
#       - the strike domain is finite
#       - the density is evaluated on a finite numerical grid
#       - the volatility surface itself is numerically fitted
#
# ============================================================================

print("\n" + "=" * 78)
print("9. DENSITY NORMALISATION CHECK")
print("=" * 78)

density_sums = np.array(
    [dbn.sum() for dbn in dbns]
)

density_errors = density_sums - 1.0

print(
    f"{'EXPIRY':>14s}"
    f"{'SUM':>16s}"
    f"{'ERROR':>16s}"
)

print("-" * 46)

for expiry_dt, density_sum, error in zip(
    expiry_dts,
    density_sums,
    density_errors,
):

    print(
        f"{str(expiry_dt):>14s}"
        f"{density_sum:16.10f}"
        f"{error:16.10f}"
    )


# Plot the density-normalisation error for each expiry.

expiry_numbers = np.arange(
    1,
    len(expiry_dts) + 1,
)

plt.figure(figsize=(10, 6))

plt.plot(
    expiry_numbers,
    density_errors,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xticks(
    expiry_numbers,
    [str(dt) for dt in expiry_dts],
    rotation=45,
)

plt.xlabel("Expiry")
plt.ylabel("Density Integral - 1")
plt.title("Risk-Neutral Density Normalisation Error")
plt.grid(True)
plt.show()


# ============================================================================
# 10. SUMMARY
# ============================================================================
#
# Summarise the main characteristics of the fitted volatility surface and the
# numerical density calculations.
#
# ============================================================================

print("\n" + "=" * 78)
print("10. SUMMARY")
print("=" * 78)

print(f"Volatility function       : {vol_function_type}")
print(f"Number of expiries        : {len(expiry_dts)}")
print(f"Number of strikes         : {len(strikes)}")
print(f"Calibration time          : {calibration_time:.6f} seconds")
print(f"Minimum density integral  : {np.min(density_sums):.10f}")
print(f"Maximum density integral  : {np.max(density_sums):.10f}")
print(
    f"Maximum normalisation err : "
    f"{np.max(np.abs(density_errors)):.10e}"
)
