# ============================================================================
# FINANCEPY EXAMPLES - Bond Portfolio
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example values a portfolio of UK government bonds using each of the
# day-count conventions supported by FinancePy.
#
# For every bond and day-count convention the example calculates:
#
#   - accrued interest
#   - yield to maturity
#
# The market clean price is held fixed. Changing the day-count convention
# changes the treatment of coupon accruals and, depending on the convention,
# can therefore change both accrued interest and the yield consistent with
# the observed clean price.
#
# This is particularly useful for illustrating why the contractual
# conventions attached to a bond are an essential part of its valuation.
# ============================================================================

import datetime as dt
import os

import pandas as pd

from financepy.products.bonds.bond import Bond
from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

set_plot_style()

# ============================================================================
# GLOBAL OUTPUT FORMAT
# ============================================================================

LINE = "=" * 110
SUBLINE = "-" * 110


# ============================================================================
# 1. LOAD GILT MARKET DATA
# ============================================================================

print("\n" + LINE)
print("1. LOAD GILT MARKET DATA")
print(LINE)

path = os.path.join(
    os.path.dirname(__file__),
    "./data/gilt_bond_prices.txt",
)

bond_dataframe = pd.read_csv(
    path,
    sep="\t",
)

bond_dataframe["mid"] = 0.5 * (bond_dataframe["bid"] + bond_dataframe["ask"])

print(f"{'Number of Bonds':<40}: " f"{len(bond_dataframe)}")

print(f"{'Data File':<40}: " f"{os.path.basename(path)}")


# ============================================================================
# 2. PORTFOLIO VALUATION SETUP
# ============================================================================
#
# The gilt portfolio is valued on 19 September 2012.
#
# All bonds are assumed to pay semi-annual coupons. The purpose of the
# example is to vary the day-count convention while keeping the remaining
# market inputs unchanged.
# ============================================================================

print("\n" + LINE)
print("2. PORTFOLIO VALUATION SETUP")
print(LINE)

settle_dt = Date(
    19,
    9,
    2012,
)

freq_type = FrequencyTypes.SEMI_ANNUAL

face = 100.0

print(f"{'Settlement Date':<40}: " f"{settle_dt}")

print(f"{'Coupon Frequency':<40}: " f"{freq_type}")

print(f"{'Face Amount':<40}: " f"{face:.2f}")


# ============================================================================
# 3. BOND PORTFOLIO ACROSS DAY-COUNT CONVENTIONS
# ============================================================================
#
# For each FinancePy day-count convention:
#
#   1. Construct every bond in the gilt portfolio.
#   2. Calculate the yield that reproduces its observed clean price.
#   3. Calculate accrued interest on a face amount of 100.
#
# DayCountTypes.ZERO is not a coupon-accrual convention and is therefore
# excluded, matching the original FinancePy test.
# ============================================================================

print("\n" + LINE)
print("3. BOND PORTFOLIO ACROSS DAY-COUNT CONVENTIONS")
print(LINE)

results = []

for dc_type in DayCountTypes:

    if dc_type == DayCountTypes.ZERO:
        continue

    print("\n" + SUBLINE)
    print(f"DAY COUNT CONVENTION: {dc_type}")
    print(SUBLINE)

    print(f"{'MATURITY':<18}" f"{'COUPON (%)':>14}" f"{'CLEAN PRICE':>16}" f"{'ACCRUED':>16}" f"{'YTM (%)':>16}")

    print(SUBLINE)

    for _, bond_row in bond_dataframe.iterrows():

        date_string = bond_row["maturity"]

        mat_date_time = dt.datetime.strptime(
            date_string,
            "%d-%b-%y",
        )

        maturity_dt = from_datetime(
            mat_date_time,
        )

        # The original FinancePy test constructs an artificial issue date
        # using the maturity day and month with the year fixed at 2000.
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

        ytm = bond.yield_to_maturity(
            settle_dt,
            clean_price,
        )

        accrued_interest = bond.accrued_interest(
            settle_dt,
            face,
        )

        print(
            f"{str(maturity_dt):<18}"
            f"{coupon * 100.0:14.6f}"
            f"{clean_price:16.6f}"
            f"{accrued_interest:16.6f}"
            f"{ytm * 100.0:16.6f}"
        )

        results.append(
            {
                "day_count": dc_type,
                "maturity": maturity_dt,
                "coupon": coupon,
                "clean_price": clean_price,
                "accrued_interest": accrued_interest,
                "ytm": ytm,
            }
        )


# ============================================================================
# 4. DAY-COUNT CONVENTION COMPARISON
# ============================================================================
#
# Select one bond from the portfolio and compare it directly across all
# day-count conventions.
#
# This isolates the effect of the convention: the bond's maturity, coupon
# and market clean price are identical in every row.
# ============================================================================

print("\n" + LINE)
print("4. DAY-COUNT CONVENTION COMPARISON")
print(LINE)

comparison_index = len(bond_dataframe) // 2

comparison_row = bond_dataframe.iloc[comparison_index]

comparison_mat_datetime = dt.datetime.strptime(
    comparison_row["maturity"],
    "%d-%b-%y",
)

comparison_maturity_dt = from_datetime(
    comparison_mat_datetime,
)

comparison_issue_dt = Date(
    comparison_maturity_dt.d,
    comparison_maturity_dt.m,
    2000,
)

comparison_coupon = comparison_row["coupon"] / 100.0

comparison_clean_price = comparison_row["mid"]


print(f"{'Selected Maturity':<40}: " f"{comparison_maturity_dt}")

print(f"{'Coupon':<40}: " f"{comparison_coupon * 100.0:12.6f}%")

print(f"{'Clean Price':<40}: " f"{comparison_clean_price:12.6f}")


print("\n" + SUBLINE)

print(f"{'DAY COUNT':<35}" f"{'ACCRUED':>18}" f"{'YTM (%)':>18}")

print(SUBLINE)


comparison_results = []

for dc_type in DayCountTypes:

    if dc_type == DayCountTypes.ZERO:
        continue

    bond = Bond(
        comparison_issue_dt,
        comparison_maturity_dt,
        comparison_coupon,
        freq_type,
        dc_type,
    )

    ytm = bond.yield_to_maturity(
        settle_dt,
        comparison_clean_price,
    )

    accrued_interest = bond.accrued_interest(
        settle_dt,
        face,
    )

    print(f"{str(dc_type):<35}" f"{accrued_interest:18.8f}" f"{ytm * 100.0:18.8f}")

    comparison_results.append(
        (
            dc_type,
            accrued_interest,
            ytm,
        )
    )


# ============================================================================
# 5. INTERPRETATION
# ============================================================================
#
# The clean market price does not change across the comparison because it is
# an observed market input.
#
# Accrued interest can change because different day-count conventions measure
# the fraction of the coupon period differently.
#
# Since:
#
#       Dirty Price = Clean Price + Accrued Interest
#
# a change in accrued interest also changes the dirty price implied by the
# same quoted clean price.
#
# Yield to maturity can consequently differ because FinancePy solves for the
# yield that reproduces the bond price under the selected contractual
# convention.
# ============================================================================

print("\n" + LINE)
print("5. INTERPRETATION")
print(LINE)

print("The market clean price is fixed for each bond.")

print("Changing the day-count convention can change accrued interest.")

print(
    "This changes the corresponding dirty price and can also change the "
    "yield to maturity consistent with the observed clean price."
)


# ============================================================================
# 6. SUMMARY
# ============================================================================

print("\n" + LINE)
print("6. SUMMARY")
print(LINE)

num_conventions = len({result["day_count"] for result in results})

print(f"{'Number of Bonds':<40}: " f"{len(bond_dataframe)}")

print(f"{'Day-Count Conventions Tested':<40}: " f"{num_conventions}")

print(f"{'Total Bond Calculations':<40}: " f"{len(results)}")

print(f"{'Settlement Date':<40}: " f"{settle_dt}")

print(f"{'Coupon Frequency':<40}: " f"{freq_type}")

print("\n" + LINE)
print("END OF BOND PORTFOLIO DEMONSTRATION")
print(LINE)
