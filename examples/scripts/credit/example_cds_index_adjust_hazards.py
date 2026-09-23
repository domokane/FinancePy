# ============================================================================
# FINANCEPY EXAMPLES - CDSIndexPortfolio
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the relationship between a CDS index and the
# constituent single-name CDS curves used to represent the index portfolio.
#
# It shows:
#
#   1. Construction of single-name issuer curves from CDS market spreads
#   2. Calculation of average single-name CDS spreads
#   3. Calculation of intrinsic index spreads
#   4. Hazard-rate adjustment of the constituent curves to index quotes
#   5. Comparison of intrinsic spreads before and after the adjustment
#
# The hazard-rate adjustment changes the constituent credit curves so that
# the intrinsic value of the portfolio is consistent with the supplied index
# market quotes.
#
# This is useful because the CDS index and its underlying single-name CDS
# contracts do not necessarily trade at exactly the same aggregate level.
# The difference is commonly associated with the index basis.
# ============================================================================

import matplotlib.pyplot as plt
import time

from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from helpers import build_ibor_curve
from helpers import load_heterogeneous_spread_curves

# ============================================================================
# OUTPUT FORMAT
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100
set_plot_style()


# ============================================================================
# 1. MARKET DATES AND INTEREST-RATE CURVE
# ============================================================================

print("\n" + LINE)
print("1. CDS INDEX PORTFOLIO")
print(LINE)

trade_dt = Date(
    1,
    8,
    2007,
)

step_in_dt = trade_dt.add_days(1)

value_dt = trade_dt

libor_curve = build_ibor_curve(trade_dt)


print(f"{'Trade Date':<30}: " f"{trade_dt}")

print(f"{'Step-In Date':<30}: " f"{step_in_dt}")

print(f"{'Value Date':<30}: " f"{value_dt}")


# ============================================================================
# 2. SINGLE-NAME CDS MATURITIES
# ============================================================================
#
# Determine standard CDS maturity dates for 3Y, 5Y, 7Y and 10Y contracts.
# ============================================================================

maturity_3yr = trade_dt.next_cds_date(36)
maturity_5yr = trade_dt.next_cds_date(60)
maturity_7yr = trade_dt.next_cds_date(84)
maturity_10yr = trade_dt.next_cds_date(120)

single_name_maturities = [
    maturity_3yr,
    maturity_5yr,
    maturity_7yr,
    maturity_10yr,
]


print("\n" + SUBLINE)

print(f"{'TENOR':>12}" f"{'MATURITY':>20}")

print(SUBLINE)

for tenor, maturity_dt in zip(
    ["3Y", "5Y", "7Y", "10Y"],
    single_name_maturities,
):
    print(f"{tenor:>12}" f"{str(maturity_dt):>20}")


# ============================================================================
# 3. LOAD SINGLE-NAME CDS MARKET DATA
# ============================================================================
#
# The input file contains CDS spreads for the constituents of CDX.NA.IG
# Series 7.
#
# Each row contains:
#
#       3Y spread
#       5Y spread
#       7Y spread
#       10Y spread
#       recovery rate
#
# CDS spreads in the file are expressed in basis points and are converted
# here to decimal form before constructing the CDS contracts.
# ============================================================================

print("\n" + LINE)
print("2. LOAD SINGLE-NAME CDS CURVES")
print(LINE)

issuer_curves = load_heterogeneous_spread_curves(value_dt, step_in_dt, libor_curve)


num_credits = len(issuer_curves)

print(f"{'Number of issuer curves':<30}: " f"{num_credits}")


# ============================================================================
# 4. AVERAGE SINGLE-NAME CDS SPREADS
# ============================================================================
#
# average_spread() calculates the average spread across the constituent
# single-name credit curves for a given maturity.
#
# This is a simple average of the individual credit spread levels and is
# distinct from the intrinsic spread of the CDS index portfolio.
# ============================================================================

print("\n" + LINE)
print("3. AVERAGE SINGLE-NAME CDS SPREADS")
print(LINE)

cds_index = CDSIndexPortfolio()


average_spreads = []

for maturity_dt in single_name_maturities:

    spread = cds_index.average_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    average_spreads.append(spread * 10000.0)


print(f"{'TENOR':>12}" f"{'AVERAGE SPREAD (bp)':>24}")

print(SUBLINE)


for tenor, spread in zip(
    ["3Y", "5Y", "7Y", "10Y"],
    average_spreads,
):

    print(f"{tenor:>12}" f"{spread:24.6f}")


# ============================================================================
# 5. INTRINSIC CDS INDEX SPREADS
# ============================================================================
#
# The intrinsic spread is calculated from the complete portfolio of
# constituent credit curves.
#
# Unlike a simple arithmetic average, the intrinsic spread reflects the
# premium and protection legs of the index portfolio.
# ============================================================================

print("\n" + LINE)
print("4. INTRINSIC CDS INDEX SPREADS")
print(LINE)


intrinsic_spreads = []

for maturity_dt in single_name_maturities:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    intrinsic_spreads.append(spread * 10000.0)


print(f"{'TENOR':>12}" f"{'AVERAGE (bp)':>20}" f"{'INTRINSIC (bp)':>20}" f"{'DIFFERENCE (bp)':>20}")

print(SUBLINE)


for tenor, average, intrinsic in zip(
    ["3Y", "5Y", "7Y", "10Y"],
    average_spreads,
    intrinsic_spreads,
):

    difference = intrinsic - average

    print(f"{tenor:>12}" f"{average:20.6f}" f"{intrinsic:20.6f}" f"{difference:20.6f}")


# ============================================================================
# 6. INDEX MARKET QUOTES
# ============================================================================
#
# These are the market index coupons and upfront payments to which the
# constituent curves will be adjusted.
#
# In this example the upfront payments are all zero. Consequently, the
# supplied index coupons represent the target index spread levels.
# ============================================================================

print("\n" + LINE)
print("5. INDEX MARKET QUOTES")
print(LINE)


index_cpns = [
    0.0020,
    0.0037,
    0.0050,
    0.0063,
]

index_upfronts = [
    0.0,
    0.0,
    0.0,
    0.0,
]

index_maturity_dts = [
    Date(20, 12, 2009),
    Date(20, 12, 2011),
    Date(20, 12, 2013),
    Date(20, 12, 2016),
]

index_recovery_rate = 0.40


print(f"{'MATURITY':>20}" f"{'COUPON (bp)':>20}" f"{'UPFRONT':>20}")

print(SUBLINE)


for maturity_dt, coupon, upfront in zip(
    index_maturity_dts,
    index_cpns,
    index_upfronts,
):

    print(f"{str(maturity_dt):>20}" f"{coupon * 10000.0:20.6f}" f"{upfront:20.6f}")


# ============================================================================
# 7. UNADJUSTED INTRINSIC SPREADS AT INDEX MATURITIES
# ============================================================================
#
# Before performing the hazard-rate adjustment, calculate the intrinsic
# spreads of the original constituent curves at exactly the same maturity
# dates as the index market quotes.
#
# This gives the starting point against which the adjustment can be judged.
# ============================================================================

print("\n" + LINE)
print("6. UNADJUSTED INTRINSIC SPREADS AT INDEX MATURITIES")
print(LINE)


unadjusted_intrinsic_spreads = []

for maturity_dt in index_maturity_dts:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    unadjusted_intrinsic_spreads.append(spread * 10000.0)


print(f"{'MATURITY':>20}" f"{'INDEX QUOTE (bp)':>22}" f"{'INTRINSIC (bp)':>22}" f"{'BASIS (bp)':>18}")

print(SUBLINE)


for maturity_dt, coupon, intrinsic in zip(
    index_maturity_dts,
    index_cpns,
    unadjusted_intrinsic_spreads,
):

    index_quote = coupon * 10000.0

    basis = intrinsic - index_quote

    print(f"{str(maturity_dt):>20}" f"{index_quote:22.6f}" f"{intrinsic:22.6f}" f"{basis:18.6f}")


# ============================================================================
# 8. HAZARD-RATE ADJUSTMENT
# ============================================================================
#
# hazard_rate_adjust_intrinsic() adjusts the constituent issuer curves so
# that the intrinsic CDS index valuation is consistent with the supplied
# index market quotes.
#
# The individual single-name curves therefore retain their role as the
# building blocks of the portfolio, but their hazard rates are adjusted to
# reconcile the portfolio with the observed index market.
# ============================================================================

print("\n" + LINE)
print("7. PERFORM CDS INDEX HAZARD-RATE ADJUSTMENT")
print(LINE)


tolerance = 1e-4

start = time.perf_counter()

adjusted_issuer_curves = cds_index.hazard_rate_adjust_intrinsic(
    value_dt,
    issuer_curves,
    index_cpns,
    index_upfronts,
    index_maturity_dts,
    index_recovery_rate,
    tolerance,
)

elapsed = time.perf_counter() - start


print(f"{'Number of issuer curves':<30}: " f"{len(adjusted_issuer_curves)}")
print(f"{'Tolerance':<30}: " f"{tolerance:.8f}")
print(f"{'Calculation time (seconds)':<30}: " f"{elapsed:.6f}")


# ============================================================================
# 9. ADJUSTED INTRINSIC SPREADS
# ============================================================================
#
# Recalculate the intrinsic index spreads using the adjusted issuer curves.
#
# If the adjustment has worked correctly, the adjusted intrinsic spreads
# should be close to the index market quotes.
# ============================================================================

print("\n" + LINE)
print("8. ADJUSTED INTRINSIC CDS INDEX SPREADS")
print(LINE)


adjusted_intrinsic_spreads = []

for maturity_dt in index_maturity_dts:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        adjusted_issuer_curves,
    )

    adjusted_intrinsic_spreads.append(spread * 10000.0)


print(f"{'MATURITY':>20}" f"{'INDEX QUOTE (bp)':>22}" f"{'ADJUSTED (bp)':>22}" f"{'ERROR (bp)':>18}")

print(SUBLINE)


for maturity_dt, coupon, adjusted in zip(
    index_maturity_dts,
    index_cpns,
    adjusted_intrinsic_spreads,
):

    index_quote = coupon * 10000.0

    error = adjusted - index_quote

    print(f"{str(maturity_dt):>20}" f"{index_quote:22.6f}" f"{adjusted:22.6f}" f"{error:18.6f}")


# ============================================================================
# 10. BEFORE AND AFTER HAZARD-RATE ADJUSTMENT
# ============================================================================
#
# Compare:
#
#       index market quote
#       original intrinsic spread
#       adjusted intrinsic spread
#
# The final column shows the residual difference between the adjusted
# intrinsic spread and the index quote.
# ============================================================================

print("\n" + LINE)
print("9. BEFORE AND AFTER HAZARD-RATE ADJUSTMENT")
print(LINE)


print(f"{'MATURITY':>20}" f"{'INDEX (bp)':>18}" f"{'BEFORE (bp)':>18}" f"{'AFTER (bp)':>18}" f"{'ERROR (bp)':>18}")

print(SUBLINE)


for (
    maturity_dt,
    coupon,
    before,
    after,
) in zip(
    index_maturity_dts,
    index_cpns,
    unadjusted_intrinsic_spreads,
    adjusted_intrinsic_spreads,
):

    index_quote = coupon * 10000.0

    error = after - index_quote

    print(f"{str(maturity_dt):>20}" f"{index_quote:18.6f}" f"{before:18.6f}" f"{after:18.6f}" f"{error:18.6f}")


# ============================================================================
# PLOT INDEX QUOTES VERSUS INTRINSIC SPREADS
# ============================================================================


index_quotes = [coupon * 10000.0 for coupon in index_cpns]

maturity_labels = [str(maturity_dt) for maturity_dt in index_maturity_dts]

plt.figure(figsize=(10, 6))

plt.plot(
    maturity_labels,
    index_quotes,
    marker="o",
    linewidth=2,
    label="Index Market Quote",
)

plt.plot(
    maturity_labels,
    unadjusted_intrinsic_spreads,
    marker="o",
    linewidth=2,
    label="Intrinsic Before Adjustment",
)

plt.plot(
    maturity_labels,
    adjusted_intrinsic_spreads,
    marker="o",
    linewidth=2,
    label="Intrinsic After Adjustment",
)

plt.xlabel("Index Maturity")
plt.ylabel("Spread (bp)")
plt.title("CDS Index Hazard-Rate Adjustment")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================================
# 11. ADJUSTMENT APPLIED TO THE INTRINSIC SPREAD
# ============================================================================
#
# This table shows how far the intrinsic portfolio spread moved as a result
# of the hazard-rate adjustment.
#
# A positive change means the adjustment increased the intrinsic spread.
# A negative change means it reduced the intrinsic spread.
# ============================================================================

print("\n" + LINE)
print("10. CHANGE IN INTRINSIC SPREAD")
print(LINE)


print(f"{'MATURITY':>20}" f"{'BEFORE (bp)':>20}" f"{'AFTER (bp)':>20}" f"{'CHANGE (bp)':>20}")

print(SUBLINE)


for maturity_dt, before, after in zip(
    index_maturity_dts,
    unadjusted_intrinsic_spreads,
    adjusted_intrinsic_spreads,
):

    change = after - before

    print(f"{str(maturity_dt):>20}" f"{before:20.6f}" f"{after:20.6f}" f"{change:20.6f}")


# ============================================================================
# 12. SUMMARY
# ============================================================================

print("\n" + LINE)
print("11. SUMMARY")
print(LINE)

print("The constituent single-name CDS spreads are first converted into " "individual issuer survival curves.")

print("Those issuer curves are then combined to calculate the intrinsic " "spread of the CDS index portfolio.")

print("The intrinsic spread need not equal the observed index market quote.")

print(
    "The hazard-rate adjustment modifies the constituent credit curves "
    "so that the intrinsic portfolio valuation is consistent with the "
    "supplied index quotes."
)

print(
    "The final comparison shows the original intrinsic spread, the "
    "adjusted intrinsic spread and the residual calibration error."
)

print("\n" + LINE)
print("END OF CDS INDEX PORTFOLIO EXAMPLE")
print(LINE)
