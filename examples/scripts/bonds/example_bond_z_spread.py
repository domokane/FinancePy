# ============================================================================
# FINANCEPY EXAMPLES - Bond Z-Spread and Asset-Swap Spread
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the calculation of:
#
#   - Z-spread
#   - Asset-swap spread
#
# for a portfolio of UK government bonds.
#
# Two base discount curves are considered:
#
#   1. A flat 1% discount curve
#   2. A market curve constructed from GBP money-market and swap benchmarks
#
# The Z-spread is the constant spread added to the discount curve that
# reproduces the observed bond price.
#
# The asset-swap spread provides an alternative measure of the bond's
# spread relative to the interest-rate curve.
#
# Results are reported in basis points and plotted against bond maturity.
# ============================================================================

from helpers import parse_date_or_tenor
import os

import matplotlib.pyplot as plt
import pandas as pd

from financepy.utils.calendar import CalendarTypes
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_vars import G_PERCENT
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.flat_discount_curve import FlatDiscountCurve
from financepy.market.curves.ibor_single_curve import IborSingleCurve
from financepy.market.curves.interpolator import InterpTypes

from financepy.products.bonds.bond import Bond
from financepy.products.rates.ibor_benchmarks_report import (
    dataframe_to_benchmarks,
)

set_plot_style()

# ============================================================================
# SUPPORTING FUNCTIONS
# ============================================================================


def calculate_bond_spreads(
    base_curve: DiscountCurve,
    curve_name: str,
):
    """Calculate Z-spreads and asset-swap spreads for the gilt portfolio."""

    # ========================================================================
    # LOAD GILT MARKET DATA
    # ========================================================================

    path = os.path.join(
        os.path.dirname(__file__),
        "./data/gilt_bond_prices.txt",
    )

    bond_dataframe = pd.read_csv(
        path,
        sep="\t",
    )

    bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

    bond_dataframe["maturity"] = pd.to_datetime(
        bond_dataframe["maturity"],
        format="%d-%b-%y",
    )

    # ========================================================================
    # BOND CONVENTIONS
    # ========================================================================

    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.ACT_ACT_ICMA

    # ========================================================================
    # CALCULATE SPREADS
    # ========================================================================
    #
    # Each bond is constructed from its maturity, coupon and standard gilt
    # conventions. The observed mid-market clean price is then used to solve
    # for both the Z-spread and asset-swap spread relative to the supplied
    # base curve.
    # ========================================================================

    for index, bond_row in bond_dataframe.iterrows():

        mat_datetime = bond_row["maturity"]

        maturity_dt = from_datetime(
            mat_datetime,
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
            accrual_type,
        )

        z_spread = bond.z_spread(
            base_curve.anchor_dt,
            clean_price,
            base_curve,
        )

        asset_swap_spread = bond.asset_swap_spread(
            base_curve.anchor_dt,
            clean_price,
            base_curve,
        )

        bond_dataframe.loc[
            index,
            "z_spread",
        ] = z_spread

        bond_dataframe.loc[
            index,
            "asset_swap_spread",
        ] = asset_swap_spread

    # ========================================================================
    # OUTPUT RESULTS
    # ========================================================================
    #
    # FinancePy returns spreads in decimal form. Multiplying by 10,000
    # converts the spreads to basis points.
    # ========================================================================

    print(f"\nCurve      : {curve_name}")
    print(f"Curve Date : {base_curve.anchor_dt}")

    print("\n" + "-" * 100)

    print(
        f"{'MATURITY':<18}"
        f"{'COUPON (%)':>14}"
        f"{'CLEAN PRICE':>16}"
        f"{'Z-SPREAD (bp)':>18}"
        f"{'ASW SPREAD (bp)':>18}"
    )

    print("-" * 100)

    for _, bond_row in bond_dataframe.iterrows():

        maturity_dt = from_datetime(
            bond_row["maturity"],
        )

        print(
            f"{str(maturity_dt):<18}"
            f"{bond_row['coupon']:14.6f}"
            f"{bond_row['mid']:16.6f}"
            f"{bond_row['z_spread'] * 10000.0:18.6f}"
            f"{bond_row['asset_swap_spread'] * 10000.0:18.6f}"
        )

    print("-" * 100)

    # ========================================================================
    # PLOT RESULTS
    # ========================================================================
    #
    # The market-data file also contains the gross redemption yield.
    #
    # Yield is already expressed in percent in the source data. Z-spread and
    # asset-swap spread are stored in decimal form, so multiplying them by
    # 100 converts them to percent for plotting on the same scale.
    # ========================================================================

    plt.figure()

    plt.plot(
        bond_dataframe["maturity"],
        bond_dataframe["gross redemption yield"],
        ".",
        label="Yield",
    )

    plt.plot(
        bond_dataframe["maturity"],
        bond_dataframe["z_spread"] * 100.0,
        ".",
        label="Z-Spread",
    )

    plt.plot(
        bond_dataframe["maturity"],
        bond_dataframe["asset_swap_spread"] * 100.0,
        ".",
        label="Asset-Swap Spread",
    )

    plt.title(f"Gilt Yield and Spread Term Structure - {curve_name}")

    plt.xlabel("Maturity")
    plt.ylabel("Percent")

    plt.legend(
        loc="best",
    )

    plt.grid()

    # ========================================================================
    # BASIC RESULT CHECKS
    # ========================================================================

    assert not bond_dataframe["z_spread"].isnull().values.any()

    assert not bond_dataframe["asset_swap_spread"].isnull().values.any()

    return bond_dataframe


# ============================================================================
# 1. Z-SPREAD USING A FLAT DISCOUNT CURVE
# ============================================================================
#
# Calculate the Z-spread and asset-swap spread for each gilt relative to a
# flat 1% discount curve.
# ============================================================================

print("\n" + "=" * 100)
print("1. Z-SPREAD USING A FLAT DISCOUNT CURVE")
print("=" * 100)

settle_dt = Date(
    19,
    9,
    2012,
)

flat_curve = FlatDiscountCurve(
    settle_dt,
    flat_zero_rate=1.0 * G_PERCENT,
)

flat_results = calculate_bond_spreads(
    flat_curve,
    "Flat 1% Curve",
)


# ============================================================================
# 2. Z-SPREAD USING A MARKET DISCOUNT CURVE
# ============================================================================
#
# Construct a GBP interest-rate curve from market benchmark instruments and
# calculate the same bond spreads relative to the resulting market curve.
# ============================================================================

print("\n" + "=" * 100)
print("2. Z-SPREAD USING A MARKET DISCOUNT CURVE")
print("=" * 100)


# ============================================================================
# 2.1 LOAD GBP CURVE MARKET DATA
# ============================================================================
#
# The benchmark file contains the deposits, FRAs and swaps used to construct
# the GBP interest-rate curve.
#
# Some start and maturity values can be tenor strings rather than explicit
# dates. The parse_date_or_tenor() helper therefore converts genuine dates
# while leaving tenor strings unchanged.
# ============================================================================

path = os.path.join(
    os.path.dirname(__file__),
    "./data/GBP_OIS_20120919.csv",
)

benchmark_dataframe = pd.read_csv(
    path,
    index_col=0,
)

benchmark_dataframe["base_date"] = benchmark_dataframe["base_date"].apply(parse_date_or_tenor)

benchmark_dataframe["start_dt"] = benchmark_dataframe["start_dt"].apply(parse_date_or_tenor)

benchmark_dataframe["maturity_dt"] = benchmark_dataframe["maturity_dt"].apply(parse_date_or_tenor)


# ============================================================================
# 2.2 CONVERT MARKET DATA TO BENCHMARK INSTRUMENTS
# ============================================================================
#
# Convert the rows in the benchmark DataFrame into FinancePy deposit, FRA
# and swap objects.
# ============================================================================

value_dt = from_datetime(
    benchmark_dataframe.loc[0, "base_date"],
)

calendar_type = CalendarTypes.LONDON

benchmarks = dataframe_to_benchmarks(
    benchmark_dataframe,
    asof_date=value_dt,
    calendar_type=calendar_type,
)

deposits = benchmarks["IborDeposit"]
fras = benchmarks["IborFRA"]
swaps = benchmarks["IborSwap"]

# Ensure the FRA instruments are ordered by maturity before curve
# construction.

fras.sort(
    key=lambda fra: fra.maturity_dt,
)


# ============================================================================
# 2.3 BUILD MARKET CURVE
# ============================================================================
#
# Bootstrap the discount curve from deposits, FRAs and swaps using linear
# interpolation of zero rates.
# ============================================================================

market_curve = IborSingleCurve(
    value_dt,
    deposits,
    fras,
    swaps,
    InterpTypes.LINEAR_ZERO_RATES,
)


# ============================================================================
# 2.4 CALCULATE BOND SPREADS
# ============================================================================
#
# Repeat the portfolio spread calculation using the bootstrapped GBP market
# curve as the reference discount curve.
# ============================================================================

market_results = calculate_bond_spreads(
    market_curve,
    "GBP Market Curve",
)


# ============================================================================
# 3. DISPLAY PLOTS
# ============================================================================
#
# Both figures have been constructed above. Calling plt.show() once at the
# end displays all plots without interrupting execution between the two
# calculations.
# ============================================================================

plt.show()
