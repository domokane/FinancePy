# ============================================================================
# FINANCEPY EXAMPLES - EquityVarianceSwap
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of an equity variance swap using
# option replication.
#
# A variance swap pays according to the difference between realised variance
# and a pre-agreed variance strike. The fair variance strike can be obtained
# from a strip of European options across strikes.
#
# This example investigates:
#
#   1. The baseline variance-swap calculation
#   2. The input volatility skew
#   3. Convergence with the number of replication options
#   4. Sensitivity to strike spacing
#   5. Sensitivity to volatility skew
#   6. Sensitivity to ATM volatility
#   7. Accuracy of the Derman skew approximation
#
# Note that variance and volatility are different quantities:
#
#       variance = volatility^2
#
# Thus a variance strike of 0.09 corresponds to a volatility strike of 30%.
#
# ============================================================================

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.market.volatility.equity_vol_curve import EquityVolCurve
from financepy.products.equity.equity_variance_swap import EquityVarianceSwap
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve

from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# HELPER FUNCTION
# ============================================================================


def vol_skew(k, atm_vol, atm_k, skew):
    """Simple linear implied-volatility skew.

    Parameters
    ----------
    k : float or ndarray
        Option strike.

    atm_vol : float
        ATM implied volatility.

    atm_k : float
        ATM strike.

    skew : float
        Change in implied volatility per unit change in strike.

    Returns
    -------
    float or ndarray
        Implied volatility at strike k.
    """

    return atm_vol + skew * (k - atm_k)


# ============================================================================
# MARKET AND CONTRACT INPUTS
# ============================================================================

start_dt = Date(20, 3, 2018)
value_dt = Date(20, 3, 2018)

tenor = "3M"

# Variance swap strikes are expressed in variance units.
#
# A 30% volatility strike therefore corresponds to:
#
#       0.30^2 = 0.09

strike = 0.30 * 0.30

vol_swap = EquityVarianceSwap(
    start_dt,
    tenor,
    strike,
)

# Equity market inputs

stock_price = 100.0
interest_rate = 0.05
dividend_yield = 0.00

discount_curve = FlatDiscountCurve(
    value_dt,
    interest_rate,
)

dividend_curve = FlatDiscountCurve(
    value_dt,
    dividend_yield,
)

# The variance swap has a three-month maturity.

t = 0.25

# Volatility-surface assumptions

atm_vol = 0.20
atm_k = 100.0

# A negative skew means that lower-strike options have higher implied
# volatility than higher-strike options.

skew = -0.02 / 5.0

strikes = np.linspace(
    50.0,
    135.0,
    18,
)

vols = vol_skew(
    strikes,
    atm_vol,
    atm_k,
    skew,
)

vol_curve = EquityVolCurve(
    strikes,
    vols,
    stock_price,
    t,
    interest_rate,
    dividend_yield,
)


# ============================================================================
# 1. BASELINE EQUITY VARIANCE SWAP
# ============================================================================

print("\n" + "=" * 78)
print("1. EQUITY VARIANCE SWAP")
print("=" * 78)

strike_spacing = 5.0

num_call_options = 10
num_put_options = 10

use_forward = False

replication_variance = vol_swap.fair_strike(
    value_dt,
    stock_price,
    dividend_curve,
    vol_curve,
    num_call_options,
    num_put_options,
    strike_spacing,
    discount_curve,
    use_forward,
)

approx_variance = vol_swap.fair_strike_approx(
    value_dt,
    stock_price,
    strikes,
    vols,
)

print(f"Contract variance strike   : {strike:.8f}")
print(f"Contract volatility strike : {np.sqrt(strike):.4%}")
print()

print(f"Replication fair variance  : {replication_variance:.8f}")
print(f"Replication fair vol       : {np.sqrt(replication_variance):.4%}")
print()

print(f"Derman approximation       : {approx_variance:.8f}")
print(f"Derman approximate vol     : {np.sqrt(approx_variance):.4%}")
print()

print(f"Approx - replication       : "
      f"{approx_variance - replication_variance:+.8f}")


# ============================================================================
# 2. PLOT THE INPUT VOLATILITY SKEW
# ============================================================================

print("\n" + "=" * 78)
print("2. INPUT VOLATILITY SKEW")
print("=" * 78)

print(f"ATM strike                 : {atm_k:.2f}")
print(f"ATM volatility             : {atm_vol:.4%}")
print(f"Volatility skew            : {skew:.6f}")

plt.figure(figsize=(10, 6))

plt.plot(
    strikes,
    100.0 * vols,
    marker="o",
    label="Implied Volatility",
)

plt.axvline(
    atm_k,
    linestyle="--",
    label="ATM Strike",
)

plt.axhline(
    100.0 * atm_vol,
    linestyle=":",
    label="ATM Volatility",
)

plt.xlabel("Strike")
plt.ylabel("Implied Volatility (%)")
plt.title("Equity Volatility Skew")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 3. NUMBER-OF-OPTIONS CONVERGENCE
# ============================================================================
#
# Variance replication theoretically involves a continuum of option strikes.
#
# Numerically we replace that continuum with a finite number of calls and
# puts. This experiment shows how the result changes as more options are
# included in the replication portfolio.
#
# ============================================================================

print("\n" + "=" * 78)
print("3. NUMBER-OF-OPTIONS CONVERGENCE")
print("=" * 78)

num_options_list = [
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
]

option_count_variances = []

print(
    f"{'OPTIONS/SIDE':>14s}"
    f"{'FAIR VARIANCE':>18s}"
    f"{'FAIR VOL':>16s}"
)

print("-" * 48)

for num_options in num_options_list:

    fair_variance = vol_swap.fair_strike(
        value_dt,
        stock_price,
        dividend_curve,
        vol_curve,
        num_options,
        num_options,
        strike_spacing,
        discount_curve,
        use_forward,
    )

    option_count_variances.append(fair_variance)

    print(
        f"{num_options:14d}"
        f"{fair_variance:18.8f}"
        f"{np.sqrt(fair_variance):16.4%}"
    )


plt.figure(figsize=(10, 6))

plt.plot(
    num_options_list,
    option_count_variances,
    marker="o",
    label="Replication",
)

plt.axhline(
    replication_variance,
    linestyle="--",
    label="Baseline",
)

plt.xlabel("Number of Calls / Puts")
plt.ylabel("Fair Variance Strike")
plt.title("Variance Swap: Number-of-Options Convergence")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. STRIKE-SPACING CONVERGENCE
# ============================================================================
#
# When changing strike spacing, also change the number of options so that
# approximately the same total strike range is covered. This separates the
# effect of discretisation from the effect of truncating the replication
# portfolio.
#
# The simple linear volatility skew used in this example becomes zero at
# K = 150, so we deliberately keep the replication range comfortably inside
# that boundary.
# ============================================================================

print("\n" + "=" * 78)
print("4. STRIKE-SPACING CONVERGENCE")
print("=" * 78)

spacing_list = [
    10.0,
    5.0,
    2.5,
    2.0,
    1.0,
]

# Keep approximately 40 strike points of coverage on each side.
replication_width = 40.0

spacing_variances = []

print(
    f"{'SPACING':>12s}"
    f"{'OPTIONS/SIDE':>16s}"
    f"{'FAIR VARIANCE':>18s}"
    f"{'FAIR VOL':>16s}"
)
print("-" * 62)

for spacing in spacing_list:

    num_options = int(replication_width / spacing)

    fair_variance = vol_swap.fair_strike(
        value_dt,
        stock_price,
        dividend_curve,
        vol_curve,
        num_options,
        num_options,
        spacing,
        discount_curve,
        use_forward,
    )

    spacing_variances.append(fair_variance)

    print(
        f"{spacing:12.2f}"
        f"{num_options:16d}"
        f"{fair_variance:18.8f}"
        f"{np.sqrt(fair_variance):16.4%}"
    )

plt.figure(figsize=(10, 6))

plt.plot(
    spacing_list,
    spacing_variances,
    marker="o",
)

plt.xlabel("Strike Spacing")
plt.ylabel("Fair Variance Strike")
plt.title("Variance Swap: Strike-Spacing Convergence")
plt.grid(True)
plt.show()


# ============================================================================
# 5. VOLATILITY-SKEW SENSITIVITY
# ============================================================================
#
# Equity index volatility surfaces generally exhibit negative skew:
#
#       low strike  -> high implied volatility
#       high strike -> low implied volatility
#
# Variance swaps are sensitive to the entire volatility smile/skew rather
# than just ATM volatility.
#
# We therefore vary the skew while holding ATM volatility fixed.
#
# ============================================================================

print("\n" + "=" * 78)
print("5. VOLATILITY-SKEW SENSITIVITY")
print("=" * 78)

# The simple linear skew is only meaningful over a limited strike range.
# With 10 replication options per side and a strike spacing of 5, skews
# more negative than about -0.004 cause the extrapolated volatility in the
# upper call wing to become non-positive. We therefore stop at -0.004.
skew_list = np.linspace(
    0.0,
    -0.004,
    9,
)

skew_replication = []
skew_approximation = []

print(
    f"{'SKEW':>12s}"
    f"{'REPLICATION':>16s}"
    f"{'APPROX':>16s}"
    f"{'DIFFERENCE':>16s}"
)

print("-" * 60)

for test_skew in skew_list:

    test_vols = vol_skew(
        strikes,
        atm_vol,
        atm_k,
        test_skew,
    )

    test_curve = EquityVolCurve(
        strikes,
        test_vols,
        stock_price,
        t,
        interest_rate,
        dividend_yield,
    )

    fair_variance = vol_swap.fair_strike(
        value_dt,
        stock_price,
        dividend_curve,
        test_curve,
        num_call_options,
        num_put_options,
        strike_spacing,
        discount_curve,
        use_forward,
    )

    approx = vol_swap.fair_strike_approx(
        value_dt,
        stock_price,
        strikes,
        test_vols,
    )

    skew_replication.append(fair_variance)
    skew_approximation.append(approx)

    print(
        f"{test_skew:12.6f}"
        f"{fair_variance:16.8f}"
        f"{approx:16.8f}"
        f"{approx - fair_variance:16.8f}"
    )


plt.figure(figsize=(10, 6))

plt.plot(
    skew_list,
    skew_replication,
    marker="o",
    label="Option Replication",
)

plt.plot(
    skew_list,
    skew_approximation,
    marker="s",
    label="Derman Approximation",
)

plt.xlabel("Volatility Skew dSigma/dK")
plt.ylabel("Fair Variance Strike")
plt.title("Variance Swap: Sensitivity to Volatility Skew")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 6. ATM-VOLATILITY SENSITIVITY
# ============================================================================
#
# Vary the overall volatility level while preserving the relative shape of
# the volatility skew.
#
# If the absolute skew were held fixed while ATM volatility were reduced,
# the simple linear volatility model could eventually generate negative
# volatilities in the upper strike wing. Instead, scale the skew with the
# ATM volatility level.
#
# The baseline parameters are:
#
#       ATM volatility = 20%
#       skew           = -0.004
#
# so each test skew is scaled in proportion to the test ATM volatility.
# ============================================================================

print("\n" + "=" * 78)
print("6. ATM-VOLATILITY SENSITIVITY")
print("=" * 78)

atm_vol_list = np.linspace(
    0.10,
    0.40,
    13,
)

atm_replication = []
atm_approximation = []

print(
    f"{'ATM VOL':>12s}"
    f"{'SKEW':>12s}"
    f"{'REPLICATION':>16s}"
    f"{'APPROX':>16s}"
    f"{'FAIR VOL':>16s}"
)

print("-" * 72)

for test_atm_vol in atm_vol_list:

    # Preserve the relative shape of the volatility skew.
    test_skew = skew * test_atm_vol / atm_vol

    test_vols = vol_skew(
        strikes,
        test_atm_vol,
        atm_k,
        test_skew,
    )

    test_curve = EquityVolCurve(
        strikes,
        test_vols,
        stock_price,
        t,
        interest_rate,
        dividend_yield,
    )

    fair_variance = vol_swap.fair_strike(
        value_dt,
        stock_price,
        dividend_curve,
        test_curve,
        num_call_options,
        num_put_options,
        strike_spacing,
        discount_curve,
        use_forward,
    )

    approx = vol_swap.fair_strike_approx(
        value_dt,
        stock_price,
        strikes,
        test_vols,
    )

    atm_replication.append(fair_variance)
    atm_approximation.append(approx)

    print(
        f"{test_atm_vol:12.4%}"
        f"{test_skew:12.6f}"
        f"{fair_variance:16.8f}"
        f"{approx:16.8f}"
        f"{np.sqrt(fair_variance):16.4%}"
    )


# ============================================================================
# 7. APPROXIMATION ERROR VERSUS SKEW
# ============================================================================
#
# Finally, look directly at the error:
#
#       approximation - replication
#
# This is useful because the approximation may work very well near a flat
# volatility surface but become progressively less accurate as the skew
# becomes stronger.
#
# ============================================================================

print("\n" + "=" * 78)
print("7. DERMAN APPROXIMATION ERROR")
print("=" * 78)

skew_replication = np.asarray(skew_replication)
skew_approximation = np.asarray(skew_approximation)

approximation_error = (
    skew_approximation - skew_replication
)

print(
    f"{'SKEW':>12s}"
    f"{'ERROR':>18s}"
)

print("-" * 30)

for test_skew, error in zip(
    skew_list,
    approximation_error,
):
    print(
        f"{test_skew:12.6f}"
        f"{error:18.8f}"
    )


plt.figure(figsize=(10, 6))

plt.plot(
    skew_list,
    approximation_error,
    marker="o",
    label="Approximation - Replication",
)

plt.axhline(
    0.0,
    linestyle="--",
    label="Zero Error",
)

plt.xlabel("Volatility Skew dSigma/dK")
plt.ylabel("Variance-Strike Error")
plt.title("Derman Approximation Error versus Volatility Skew")
plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 8. SUMMARY
# ============================================================================

print("\n" + "=" * 78)
print("8. SUMMARY")
print("=" * 78)

print(f"ATM volatility             : {atm_vol:.4%}")
print(f"Volatility skew            : {skew:.6f}")
print()
print(f"Contract variance strike   : {strike:.8f}")
print(f"Contract volatility strike : {np.sqrt(strike):.4%}")
print()
print(f"Replication fair variance  : {replication_variance:.8f}")
print(f"Replication fair vol       : {np.sqrt(replication_variance):.4%}")
print()
print(f"Derman approximation       : {approx_variance:.8f}")
print(f"Derman approximate vol     : {np.sqrt(approx_variance):.4%}")
print()
print(
    f"Approximation error        : "
    f"{approx_variance - replication_variance:+.8f}"
)
