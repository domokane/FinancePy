# ============================================================================
# FINANCEPY EXAMPLES - CDSIndexPortfolio
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates spread adjustment of the constituent issuer
# curves of a CDS index portfolio.
#
# It shows:
#
#   1. Construction of the interest-rate curve
#   2. Construction of single-name CDS issuer curves
#   3. Average constituent CDS spreads
#   4. Intrinsic CDS index spreads
#   5. Spread adjustment to observed index market quotes
#   6. Comparison of intrinsic spreads before and after adjustment
#
# The spread adjustment modifies the constituent issuer curves so that the
# intrinsic index spreads are consistent with the supplied index market
# quotes.
# ============================================================================

import time

import matplotlib.pyplot as plt

from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

from helpers import build_ibor_curve
from helpers import load_heterogeneous_issuer_curves

# ============================================================================
# 1. CDS INDEX ADJUST SPREADS
# ============================================================================

print("\n" + "=" * 100)
print("1. CDS INDEX ADJUST SPREADS")
print("=" * 100)
set_plot_style()


# ============================================================================
# 1.1 MARKET DATES
# ============================================================================

trade_dt = Date(1, 8, 2007)
step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt

libor_curve = build_ibor_curve(
    trade_dt,
)

print(f"{'Trade Date':<30}: {trade_dt}")
print(f"{'Step-In Date':<30}: {step_in_dt}")
print(f"{'Value Date':<30}: {value_dt}")


# ============================================================================
# 1.2 STANDARD SINGLE-NAME CDS MATURITIES
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


# ============================================================================
# 1.3 LOAD SINGLE-NAME CDS MARKET DATA
# ============================================================================

issuer_curves = load_heterogeneous_issuer_curves(value_dt, step_in_dt, libor_curve)

print(f"{'Number of issuer curves':<30}: {len(issuer_curves)}")


# ============================================================================
# 2. AVERAGE SINGLE-NAME CDS SPREADS
# ============================================================================
#
# The average spread is the average spread across the constituent issuer
# curves. It is not necessarily equal to the intrinsic spread of the index.
# ============================================================================

print("\n" + "=" * 100)
print("2. AVERAGE SINGLE-NAME CDS SPREADS")
print("=" * 100)

cds_index = CDSIndexPortfolio()

average_spreads = []

for maturity_dt in single_name_maturities:

    spread = cds_index.average_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    average_spreads.append(
        spread * 10000.0,
    )


print(f"{'TENOR':>12}" f"{'MATURITY':>20}" f"{'AVERAGE SPREAD (bp)':>24}")

print("-" * 56)

for tenor, maturity_dt, spread in zip(
    ["3Y", "5Y", "7Y", "10Y"],
    single_name_maturities,
    average_spreads,
):

    print(f"{tenor:>12}" f"{str(maturity_dt):>20}" f"{spread:24.6f}")


# ============================================================================
# 3. INTRINSIC CDS INDEX SPREADS
# ============================================================================
#
# The intrinsic spread is calculated from the complete portfolio of
# constituent issuer curves.
# ============================================================================

print("\n" + "=" * 100)
print("3. INTRINSIC CDS INDEX SPREADS")
print("=" * 100)

intrinsic_spreads = []

for maturity_dt in single_name_maturities:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    intrinsic_spreads.append(
        spread * 10000.0,
    )


print(f"{'TENOR':>12}" f"{'AVERAGE (bp)':>20}" f"{'INTRINSIC (bp)':>20}" f"{'DIFFERENCE (bp)':>20}")

print("-" * 72)

for tenor, average, intrinsic in zip(
    ["3Y", "5Y", "7Y", "10Y"],
    average_spreads,
    intrinsic_spreads,
):

    difference = intrinsic - average

    print(f"{tenor:>12}" f"{average:20.6f}" f"{intrinsic:20.6f}" f"{difference:20.6f}")


# ============================================================================
# 4. INDEX MARKET QUOTES
# ============================================================================

print("\n" + "=" * 100)
print("4. INDEX MARKET QUOTES")
print("=" * 100)

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

index_recovery = 0.40

index_quotes = [cpn * 10000.0 for cpn in index_cpns]


print(f"{'MATURITY':>20}" f"{'COUPON (bp)':>20}" f"{'UPFRONT':>20}")

print("-" * 60)

for maturity_dt, coupon, upfront in zip(
    index_maturity_dts,
    index_cpns,
    index_upfronts,
):

    print(f"{str(maturity_dt):>20}" f"{coupon * 10000.0:20.6f}" f"{upfront:20.6f}")


# ============================================================================
# 5. INTRINSIC SPREADS BEFORE ADJUSTMENT
# ============================================================================
#
# Calculate the intrinsic portfolio spreads at exactly the same maturity
# dates as the index market quotes.
#
# This is important because the standard single-name CDS maturities above
# are not necessarily identical to the index maturity dates.
# ============================================================================

print("\n" + "=" * 100)
print("5. INTRINSIC SPREADS BEFORE ADJUSTMENT")
print("=" * 100)

intrinsic_before = []

for maturity_dt in index_maturity_dts:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        issuer_curves,
    )

    intrinsic_before.append(
        spread * 10000.0,
    )


print(f"{'MATURITY':>20}" f"{'INDEX (bp)':>20}" f"{'INTRINSIC (bp)':>20}" f"{'BASIS (bp)':>20}")

print("-" * 80)

for maturity_dt, index_quote, intrinsic in zip(
    index_maturity_dts,
    index_quotes,
    intrinsic_before,
):

    basis = intrinsic - index_quote

    print(f"{str(maturity_dt):>20}" f"{index_quote:20.6f}" f"{intrinsic:20.6f}" f"{basis:20.6f}")


# ============================================================================
# 6. PERFORM SPREAD ADJUSTMENT
# ============================================================================
#
# spread_adjust_intrinsic() adjusts the constituent issuer curves so that
# the intrinsic portfolio spreads reproduce the supplied CDS index market
# quotes.
# ============================================================================

print("\n" + "=" * 100)
print("6. PERFORM CDS INDEX SPREAD ADJUSTMENT")
print("=" * 100)

tolerance = 1e-4

start = time.perf_counter()

index_portfolio = CDSIndexPortfolio()

adjusted_issuer_curves = index_portfolio.spread_adjust_intrinsic(
    value_dt,
    issuer_curves,
    index_cpns,
    index_upfronts,
    index_maturity_dts,
    index_recovery,
    tolerance,
)

elapsed = time.perf_counter() - start

print(f"{'Tolerance':<35}: {tolerance:.8f}")
print(f"{'Number of issuer curves':<35}: {len(adjusted_issuer_curves)}")
print(f"{'Calculation time (seconds)':<35}: {elapsed:.6f}")


# ============================================================================
# 7. INTRINSIC SPREADS AFTER ADJUSTMENT
# ============================================================================

print("\n" + "=" * 100)
print("7. INTRINSIC SPREADS AFTER ADJUSTMENT")
print("=" * 100)

intrinsic_after = []

for maturity_dt in index_maturity_dts:

    spread = cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        maturity_dt,
        adjusted_issuer_curves,
    )

    intrinsic_after.append(
        spread * 10000.0,
    )


print(f"{'MATURITY':>20}" f"{'INDEX (bp)':>20}" f"{'ADJUSTED (bp)':>20}" f"{'ERROR (bp)':>20}")

print("-" * 80)

for maturity_dt, index_quote, adjusted in zip(
    index_maturity_dts,
    index_quotes,
    intrinsic_after,
):

    error = adjusted - index_quote

    print(f"{str(maturity_dt):>20}" f"{index_quote:20.6f}" f"{adjusted:20.6f}" f"{error:20.6f}")


# ============================================================================
# 8. BEFORE AND AFTER SPREAD ADJUSTMENT
# ============================================================================

print("\n" + "=" * 100)
print("8. BEFORE AND AFTER SPREAD ADJUSTMENT")
print("=" * 100)

print(
    f"{'MATURITY':>20}"
    f"{'INDEX (bp)':>18}"
    f"{'BEFORE (bp)':>18}"
    f"{'AFTER (bp)':>18}"
    f"{'CHANGE (bp)':>18}"
    f"{'ERROR (bp)':>18}"
)

print("-" * 110)

for (
    maturity_dt,
    index_quote,
    before,
    after,
) in zip(
    index_maturity_dts,
    index_quotes,
    intrinsic_before,
    intrinsic_after,
):

    change = after - before
    error = after - index_quote

    print(
        f"{str(maturity_dt):>20}"
        f"{index_quote:18.6f}"
        f"{before:18.6f}"
        f"{after:18.6f}"
        f"{change:18.6f}"
        f"{error:18.6f}"
    )


# ============================================================================
# 9. PLOT CDS INDEX SPREAD ADJUSTMENT
# ============================================================================
#
# The chart compares:
#
#   - observed CDS index market quotes
#   - intrinsic spreads before adjustment
#   - intrinsic spreads after adjustment
#
# If the calibration is successful, the adjusted intrinsic spread should
# lie close to the corresponding index market quote.
# ============================================================================

print("\n" + "=" * 100)
print("9. PLOT CDS INDEX SPREAD ADJUSTMENT")
print("=" * 100)

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
    intrinsic_before,
    marker="o",
    linewidth=2,
    label="Intrinsic Before Adjustment",
)

plt.plot(
    maturity_labels,
    intrinsic_after,
    marker="o",
    linewidth=2,
    label="Intrinsic After Adjustment",
)

plt.xlabel("Index Maturity")
plt.ylabel("Spread (bp)")
plt.title("CDS Index Spread Adjustment")

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ============================================================================
# 10. SUMMARY
# ============================================================================

print("\n" + "=" * 100)
print("10. SUMMARY")
print("=" * 100)

print("The constituent single-name CDS market spreads are first used to " "construct individual issuer credit curves.")

print("Those issuer curves are combined to calculate the intrinsic spread " "of the CDS index portfolio.")

print("The intrinsic portfolio spread can differ from the observed CDS " "index market quote.")

print(
    "The spread adjustment modifies the constituent issuer curves so that "
    "the intrinsic portfolio spreads are consistent with the supplied "
    "index market quotes."
)

print(
    "The final table and plot compare the index quote with the intrinsic "
    "spread before and after the spread adjustment."
)

print("\n" + "=" * 100)
print("END OF CDS INDEX SPREAD ADJUSTMENT EXAMPLE")
print("=" * 100)
