# ============================================================================
# FINANCEPY EXAMPLES - BondParametricDiscountCurve
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates how a parametric discount curve can be fitted
# to a cross-section of government bond prices.
#
# Several curve-fitting methods are compared:
#
#   1. Cubic polynomial
#   2. Quintic polynomial
#   3. Nelson-Siegel
#   4. Nelson-Siegel-Svensson
#   5. B-spline
#
# For each fitted curve, the example:
#
#   - fits the curve to observed gilt bond prices
#   - plots the fitted bond yield curve
#   - calculates discount factors
#   - calculates continuously compounded zero rates
#   - calculates instantaneous forward rates
#   - reports RMS yield and price fitting errors
#
# This makes it possible to compare both the shape of the fitted curves and
# their ability to reproduce the observed bond market data.
# ============================================================================

import datetime as dt
import os

import numpy as np
import pandas as pd

from financepy.market.curves import BondParametricDiscountCurve
from financepy.market.curves.curve_fits import CurveFitTypes
from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# 1. BOND PARAMETRIC DISCOUNT CURVE
# ============================================================================

print("\n" + "=" * 100)
print("1. BOND PARAMETRIC DISCOUNT CURVE")
print("=" * 100)


# ============================================================================
# 1.1 LOAD GILT MARKET DATA
# ============================================================================
#
# Read bid and ask prices for the gilt bonds and use the mid-market price
# as the clean price supplied to the curve-fitting procedure.
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
# 1.2 BUILD BOND INSTRUMENTS
# ============================================================================
#
# Construct a Bond object for every bond in the market-data file.
#
# The example assumes semi-annual coupons and ACT/ACT ICMA accrual.
# ============================================================================

freq_type = FrequencyTypes.SEMI_ANNUAL
dc_type = DayCountTypes.ACT_ACT_ICMA

settle_dt = Date(19, 9, 2012)

bonds = []
clean_prices = []

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

    bonds.append(bond)
    clean_prices.append(clean_price)


# ============================================================================
# 1.3 CURVE-FITTING METHODS
# ============================================================================
#
# Fit the same bond data using several alternative functional forms.
#
# Comparing the resulting curves illustrates how the choice of fitting
# method affects the term structure implied by the same market prices.
# ============================================================================

curve_fitters = [
    CurveFitTypes.CUBIC_POLYNOMIAL,
    CurveFitTypes.QUINTIC_POLYNOMIAL,
    CurveFitTypes.NELSON_SIEGEL,
    CurveFitTypes.NELSON_SIEGEL_SVENSSON,
    CurveFitTypes.BSPLINE,
]


# ============================================================================
# 1.4 CURVE EVALUATION GRID
# ============================================================================
#
# Evaluate each fitted curve from approximately zero to 50 years.
#
# A small positive value is used instead of exactly zero.
# ============================================================================

times = np.linspace(
    1.0e-6,
    50.0,
    100,
)


# ============================================================================
# 1.5 FIT AND COMPARE CURVES
# ============================================================================

for curve_fitter in curve_fitters:

    name = curve_fitter.name

    print("\n" + "=" * 100)
    print(f"CURVE FITTER: {name}")
    print("=" * 100)

    # ========================================================================
    # FIT PARAMETRIC CURVE
    # ========================================================================

    fitted_curve = BondParametricDiscountCurve(
        settle_dt,
        bonds,
        clean_prices,
        curve_fitter,
    )

    # ========================================================================
    # PLOT BOND YIELD FIT
    # ========================================================================
    #
    # Plot the observed bond yields together with the fitted parametric
    # curve. The plot is generated for every curve-fitting method.
    # ========================================================================

    fitted_curve.plot_bond_yield_fit(
        name + " Bond Yield Fit",
    )

    # ========================================================================
    # CALCULATE TERM-STRUCTURE QUANTITIES
    # ========================================================================
    #
    # Evaluate:
    #
    #   df_t()             -> discount factor
    #   zero_rate_cc_t()   -> continuously compounded zero rate
    #   fwd_rate_inst_t()  -> instantaneous forward rate
    # ========================================================================

    dfs = fitted_curve.df_t(
        times,
    )

    zeros = fitted_curve.zero_rate_cc_t(
        times,
    )

    fwds = fitted_curve.fwd_rate_inst_t(
        times,
    )

    # ========================================================================
    # TERM-STRUCTURE OUTPUT
    # ========================================================================

    print("\nTERM STRUCTURE")
    print("-" * 78)

    print(f"{'TIME':>12}" f"{'DISCOUNT FACTOR':>20}" f"{'ZERO RATE':>20}" f"{'FORWARD RATE':>20}")

    print("-" * 78)

    for i, t in enumerate(times):

        print(f"{t:12.6f}" f"{dfs[i]:20.10f}" f"{zeros[i]:20.10f}" f"{fwds[i]:20.10f}")

    # ========================================================================
    # FITTING ERRORS
    # ========================================================================
    #
    # RMS yield error measures the quality of the fit in yield space.
    #
    # RMS price error measures the quality of the fit in bond-price space.
    # ========================================================================

    rms_yield_err = fitted_curve.rms_yield_error()
    rms_price_err = fitted_curve.rms_price_error()

    print("\nFIT ERRORS")
    print("-" * 58)

    print(f"{'ERROR TYPE':<30}" f"{'VALUE':>20}")

    print("-" * 58)

    print(f"{'RMS Yield Error':<30}" f"{rms_yield_err:20.10f}")

    print(f"{'RMS Price Error':<30}" f"{rms_price_err:20.10f}")

    print("-" * 58)
