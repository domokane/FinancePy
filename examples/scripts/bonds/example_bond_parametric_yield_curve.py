# ============================================================================
# FINANCEPY EXAMPLES - BondParametricYieldCurve
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates how parametric yield curves can be fitted to
# yields obtained from a cross-section of government bond prices.
#
# The procedure is:
#
#   1. Read observed gilt bond prices.
#   2. Calculate each bond's yield to maturity.
#   3. Fit several parametric yield-curve models to those yields.
#   4. Plot each fitted yield curve.
#   5. Report the fitted model parameters.
#   6. Interpolate a yield from the fitted B-spline curve.
#
# The curve-fitting methods considered are:
#
#   - Cubic polynomial
#   - Quintic polynomial
#   - Nelson-Siegel
#   - Nelson-Siegel-Svensson
#   - B-spline
# ============================================================================

import datetime as dt
import os

import pandas as pd

from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves import BondParametricYieldCurve
from financepy.market.curves import CurveFitTypes

from financepy.products.bonds.bond import Bond

set_plot_style()

# ============================================================================
# 1. BOND PARAMETRIC YIELD CURVE
# ============================================================================

print("\n" + "=" * 100)
print("1. BOND PARAMETRIC YIELD CURVE")
print("=" * 100)


# ============================================================================
# 1.1 LOAD GILT MARKET DATA
# ============================================================================
#
# Read bid and ask prices and use the mid-market price as the clean price
# for each bond.
# ============================================================================

path = os.path.join(
    os.path.dirname(__file__),
    "./data/gilt_bond_prices.txt",
)

bond_dataframe = pd.read_csv(
    path,
    sep="\t",
)

bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])


# ============================================================================
# 1.2 BUILD BONDS AND CALCULATE MARKET YIELDS
# ============================================================================
#
# Construct a Bond object for every bond in the data set and solve for the
# yield to maturity that reproduces its observed clean price.
#
# These market-implied yields are the observations used in the subsequent
# parametric yield-curve fits.
# ============================================================================

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

settle_dt = Date(19, 9, 2012)

bonds = []
ylds = []

print("\n" + "-" * 82)

print(f"{'MATURITY':<18}" f"{'COUPON (%)':>14}" f"{'CLEAN PRICE':>18}" f"{'YIELD (%)':>16}")

print("-" * 82)

for _, bond_row in bond_dataframe.iterrows():

    date_string = bond_row["maturity"]

    mat_date_time = dt.datetime.strptime(
        date_string,
        "%d-%b-%y",
    )

    maturity_dt = from_datetime(
        mat_date_time,
    )

    issue_dt = Date(
        maturity_dt.d,
        maturity_dt.m,
        2000,
    )

    coupon = bond_row["coupon"] / 100.0
    clean_price = bond_row["mid"]

    bond = Bond(
        issue_dt,
        maturity_dt,
        coupon,
        freq_type,
        dc_type,
    )

    yld = bond.yield_to_maturity(
        settle_dt,
        clean_price,
    )

    bonds.append(bond)
    ylds.append(yld)

    print(f"{str(maturity_dt):<18}" f"{coupon * 100.0:14.6f}" f"{clean_price:18.6f}" f"{yld * 100.0:16.6f}")

print("-" * 82)


# ============================================================================
# 1.3 CUBIC POLYNOMIAL
# ============================================================================

fit_type = CurveFitTypes.CUBIC_POLYNOMIAL

fitted_curve1 = BondParametricYieldCurve(
    settle_dt,
    bonds,
    ylds,
    fit_type,
)

fitted_curve1.plot(
    "GBP Yield Curve - Cubic Polynomial",
)


# ============================================================================
# 1.4 QUINTIC POLYNOMIAL
# ============================================================================

fit_type = CurveFitTypes.QUINTIC_POLYNOMIAL

fitted_curve2 = BondParametricYieldCurve(
    settle_dt,
    bonds,
    ylds,
    fit_type,
)

fitted_curve2.plot(
    "GBP Yield Curve - Quintic Polynomial",
)


# ============================================================================
# 1.5 NELSON-SIEGEL
# ============================================================================

fit_type = CurveFitTypes.NELSON_SIEGEL

fitted_curve3 = BondParametricYieldCurve(
    settle_dt,
    bonds,
    ylds,
    fit_type,
)

fitted_curve3.plot(
    "GBP Yield Curve - Nelson-Siegel",
)


# ============================================================================
# 1.6 NELSON-SIEGEL-SVENSSON
# ============================================================================

fit_type = CurveFitTypes.NELSON_SIEGEL_SVENSSON

fitted_curve4 = BondParametricYieldCurve(
    settle_dt,
    bonds,
    ylds,
    fit_type,
)

fitted_curve4.plot(
    "GBP Yield Curve - Nelson-Siegel-Svensson",
)


# ============================================================================
# 1.7 B-SPLINE
# ============================================================================

fit_type = CurveFitTypes.BSPLINE

fitted_curve5 = BondParametricYieldCurve(
    settle_dt,
    bonds,
    ylds,
    fit_type,
)

fitted_curve5.plot(
    "GBP Yield Curve - B-Spline",
)


# ============================================================================
# 2. FITTED CURVE PARAMETERS
# ============================================================================
#
# Report the parameters estimated by each parametric curve fit.
#
# Polynomial curves are represented by their fitted coefficient vectors.
#
# Nelson-Siegel and Nelson-Siegel-Svensson expose their fitted beta and tau
# parameters directly.
# ============================================================================

print("\n" + "=" * 100)
print("2. FITTED CURVE PARAMETERS")
print("=" * 100)


# ============================================================================
# 2.1 CUBIC POLYNOMIAL
# ============================================================================

print("\nCUBIC POLYNOMIAL")
print("-" * 54)

for i, coeff in enumerate(fitted_curve1.curve_fit.coeffs):
    print(f"{'Coefficient ' + str(i):<30}" f"{coeff:20.10f}")


# ============================================================================
# 2.2 QUINTIC POLYNOMIAL
# ============================================================================

print("\nQUINTIC POLYNOMIAL")
print("-" * 54)

for i, coeff in enumerate(fitted_curve2.curve_fit.coeffs):
    print(f"{'Coefficient ' + str(i):<30}" f"{coeff:20.10f}")


# ============================================================================
# 2.3 NELSON-SIEGEL
# ============================================================================

print("\nNELSON-SIEGEL")
print("-" * 54)

print(f"{'Beta 1':<30}" f"{fitted_curve3.curve_fit.beta_1:20.10f}")

print(f"{'Beta 2':<30}" f"{fitted_curve3.curve_fit.beta_2:20.10f}")

print(f"{'Beta 3':<30}" f"{fitted_curve3.curve_fit.beta_3:20.10f}")

print(f"{'Tau':<30}" f"{fitted_curve3.curve_fit.tau:20.10f}")


# ============================================================================
# 2.4 NELSON-SIEGEL-SVENSSON
# ============================================================================

print("\nNELSON-SIEGEL-SVENSSON")
print("-" * 54)

print(f"{'Beta 1':<30}" f"{fitted_curve4.curve_fit.beta_1:20.10f}")

print(f"{'Beta 2':<30}" f"{fitted_curve4.curve_fit.beta_2:20.10f}")

print(f"{'Beta 3':<30}" f"{fitted_curve4.curve_fit.beta_3:20.10f}")

print(f"{'Beta 4':<30}" f"{fitted_curve4.curve_fit.beta_4:20.10f}")

print(f"{'Tau 1':<30}" f"{fitted_curve4.curve_fit.tau_1:20.10f}")

print(f"{'Tau 2':<30}" f"{fitted_curve4.curve_fit.tau_2:20.10f}")


# ============================================================================
# 3. INTERPOLATED YIELD
# ============================================================================
#
# Interpolate the yield for a selected maturity using the fitted B-spline
# yield curve.
# ============================================================================

print("\n" + "=" * 100)
print("3. INTERPOLATED YIELD - B-SPLINE")
print("=" * 100)

maturity_dt = Date(
    19,
    9,
    2030,
)

interp_yield = fitted_curve5.interp_yield(
    maturity_dt,
)

print(f"{'Maturity Date':<30}" f"{str(maturity_dt):>20}")

print(f"{'Interpolated Yield':<30}" f"{interp_yield:20.10f}")

print(f"{'Interpolated Yield (%)':<30}" f"{interp_yield * 100.0:20.6f}")
