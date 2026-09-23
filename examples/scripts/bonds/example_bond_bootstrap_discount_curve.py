# ============================================================================
# FINANCEPY EXAMPLES - BondBootstrapDiscountCurve
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates how to bootstrap a discount curve from a set of
# observed government bond prices.
#
# The example uses UK gilt market data and:
#
#   1. Loads bid and ask prices and calculates mid-market prices
#   2. Applies UK government bond market conventions
#   3. Constructs Bond objects for the observed gilts
#   4. Bootstraps a discount curve from the bond prices
#   5. Repeats the bootstrap for each available interpolation method
#   6. Reports zero rates and discount factors at bond maturities
#   7. Compares the resulting zero-rate and forward-rate curves
#
# Bootstrapping determines discount factors that reproduce the observed
# market prices of the calibration instruments. The interpolation method
# determines how the curve behaves between the calibrated maturity points.
# ============================================================================

import datetime as dt
import os
import time

import matplotlib.pyplot as plt
import pandas as pd

from financepy.utils.date import Date
from financepy.utils.date import from_datetime
from financepy.utils.date_format import DateFormatTypes
from financepy.utils.date_format import set_date_format
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves import BondBootstrapDiscountCurve
from financepy.market.curves.interpolator import InterpTypes

from financepy.products.bonds import BondMarkets
from financepy.products.bonds import get_bond_market_conventions
from financepy.products.bonds.bond import Bond

set_date_format(DateFormatTypes.UK_LONG)
set_plot_style()


# ============================================================================
# 1. LOAD UK GILT MARKET DATA
# ============================================================================
#
# Read the observed gilt bid and ask prices used to construct the curve.
# The mid-market clean price is used as the calibration price for each bond.
# ============================================================================

print("\n" + "=" * 78)
print("1. BOND BOOTSTRAP DISCOUNT CURVE")
print("=" * 78)

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
# 2. UK GOVERNMENT BOND CONVENTIONS
# ============================================================================
#
# Obtain the standard market conventions for UK government bonds. These
# determine the accrual day-count convention, coupon frequency, settlement
# lag, ex-dividend period and calendar used when constructing each gilt.
# ============================================================================

(
    acc_dc_type,
    freq_type,
    spot_days,
    ex_div_days,
    cal,
) = get_bond_market_conventions(BondMarkets.UNITED_KINGDOM)

today = Date(18, 9, 2012)
settle_dt = today.add_weekdays(spot_days)


# ============================================================================
# 3. CONSTRUCT THE GILT CALIBRATION PORTFOLIO
# ============================================================================
#
# Convert each row of market data into a FinancePy Bond object and retain
# its observed mid-market clean price.
#
# The resulting bond portfolio provides the calibration instruments used
# by BondBootstrapDiscountCurve.
# ============================================================================

bonds = []
clean_prices = []

for _, bond_row in bond_dataframe.iterrows():

    date_string = bond_row["maturity"]
    mat_date_time = dt.datetime.strptime(
        date_string,
        "%d-%b-%y",
    )

    maturity_dt = from_datetime(mat_date_time)

    # The source data does not provide issue dates. The original example
    # assigns an issue date in 2000 while preserving the maturity day and
    # month.
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
        acc_dc_type,
        ex_div_days=ex_div_days,
        cal_type=cal,
    )

    bonds.append(bond)
    clean_prices.append(clean_price)


# ============================================================================
# 4. BOOTSTRAP USING EACH INTERPOLATION METHOD
# ============================================================================
#
# Construct a separate bond discount curve for every interpolation method
# available through InterpTypes.
#
# check_refit_flag=True checks that the bootstrapped curve can reproduce
# the prices of the bonds used in the calibration.
#
# The elapsed bootstrap time is also reported to allow the computational
# cost of the interpolation methods to be compared.
# ============================================================================

for interp_type in InterpTypes:

    start = time.perf_counter()

    bond_curve = BondBootstrapDiscountCurve(
        settle_dt,
        bonds,
        clean_prices,
        interp_type,
        check_refit_flag=True,
    )

    elapsed = time.perf_counter() - start

    print("\n" + "-" * 78)
    print("INTERPOLATION METHOD:", interp_type)
    print("BOOTSTRAP TIME:", elapsed)
    print("-" * 78)

    print(
        "DATE",
        "ZERO RATE",
        "DISCOUNT FACTOR",
    )

    # Evaluate the calibrated curve at each observed bond maturity.
    for _, bond_row in bond_dataframe.iterrows():

        date_string = bond_row["maturity"]
        mat_date_time = dt.datetime.strptime(
            date_string,
            "%d-%b-%y",
        )

        maturity_dt = from_datetime(mat_date_time)

        zero_rate = bond_curve.zero_rate(maturity_dt)
        df = bond_curve.df(maturity_dt)

        print(
            maturity_dt,
            zero_rate * 100,
            df,
        )

    # Plot the term structures implied by the current interpolation method.
    #
    # The zero-rate curve shows the continuously constructed term structure
    # of spot rates, while the forward-rate curve highlights the local shape
    # implied between maturity points.
    bond_curve.plot_zero_rates("BOND EXACT CURVE - " + str(interp_type))

    bond_curve.plot_fwd_rates("BOND EXACT CURVE - " + str(interp_type))

    plt.show()
