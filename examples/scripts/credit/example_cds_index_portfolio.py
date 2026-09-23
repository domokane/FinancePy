# ============================================================================
# FINANCEPY EXAMPLES - CDSIndexPortfolio
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

from helpers import load_heterogeneous_spread_curves
from helpers import build_ibor_curve
import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio

set_plot_style()


# ============================================================================
# 1. CDS INDEX PORTFOLIO
# ============================================================================
# What this section demonstrates:
# The loop varies dates, parameters, instruments or conventions so their effect can be compared rather than relying on one isolated result.

print("\n" + "=" * 78)
print("1. CDS INDEX PORTFOLIO")
print("=" * 78)

trade_dt = Date(1, 8, 2007)
step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt


libor_curve = build_ibor_curve(trade_dt)

maturity_3yr = trade_dt.next_cds_date(36)
maturity_5yr = trade_dt.next_cds_date(60)
maturity_7yr = trade_dt.next_cds_date(84)
maturity_10yr = trade_dt.next_cds_date(120)

issuer_curves = load_heterogeneous_spread_curves(value_dt, step_in_dt, libor_curve)

print(f"{'Number of Issuers':<30}: {len(issuer_curves)}")

# Now determine the average spread of the index

cds_index = CDSIndexPortfolio()

avg_spd_3yr = cds_index.average_spread(value_dt, step_in_dt, maturity_3yr, issuer_curves) * 10000.0

avg_spd_5yr = cds_index.average_spread(value_dt, step_in_dt, maturity_5yr, issuer_curves) * 10000.0

avg_spd_7yr = cds_index.average_spread(value_dt, step_in_dt, maturity_7yr, issuer_curves) * 10000.0

avg_spd_10yr = cds_index.average_spread(value_dt, step_in_dt, maturity_10yr, issuer_curves) * 10000.0

print("LABEL", "VALUE")
print("AVERAGE SPD 3Y", avg_spd_3yr)
print("AVERAGE SPD 5Y", avg_spd_5yr)
print("AVERAGE SPD 7Y", avg_spd_7yr)
print("AVERAGE SPD 10Y", avg_spd_10yr)

# Now determine the intrinsic spread of the index to the same maturity
# dates. As the single name CDS contracts

cds_index = CDSIndexPortfolio()

intrinsic_spd_3yr = cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_3yr, issuer_curves) * 10000.0

intrinsic_spd_5yr = cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_5yr, issuer_curves) * 10000.0

intrinsic_spd_7yr = cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_7yr, issuer_curves) * 10000.0

intrinsic_spd_10yr = cds_index.intrinsic_spread(value_dt, step_in_dt, maturity_10yr, issuer_curves) * 10000.0

print("LABEL", "VALUE")
print("INTRINSIC SPD 3Y", intrinsic_spd_3yr)
print("INTRINSIC SPD 5Y", intrinsic_spd_5yr)
print("INTRINSIC SPD 7Y", intrinsic_spd_7yr)
print("INTRINSIC SPD 10Y", intrinsic_spd_10yr)

# ============================================================================
# 2. COMPARE AVERAGE AND INTRINSIC INDEX SPREADS
# ============================================================================
#
# The average spread is the simple average of the constituent CDS spreads.
#
# The intrinsic spread is obtained from the aggregate protection and premium
# legs of the portfolio. It therefore reflects the valuation mechanics of the
# CDS index portfolio rather than simply averaging constituent spreads.
#
# Comparing the two across maturity illustrates the difference between these
# two measures of index credit spread.
# ============================================================================

print("\n" + "=" * 78)
print("2. AVERAGE VERSUS INTRINSIC INDEX SPREAD")
print("=" * 78)

tenors = np.array(
    [
        3,
        5,
        7,
        10,
    ]
)

average_spreads = np.array(
    [
        avg_spd_3yr,
        avg_spd_5yr,
        avg_spd_7yr,
        avg_spd_10yr,
    ]
)

intrinsic_spreads = np.array(
    [
        intrinsic_spd_3yr,
        intrinsic_spd_5yr,
        intrinsic_spd_7yr,
        intrinsic_spd_10yr,
    ]
)

print(f"{'TENOR':>10}" f"{'AVERAGE (bp)':>20}" f"{'INTRINSIC (bp)':>20}" f"{'DIFFERENCE (bp)':>20}")

print("-" * 70)

for tenor, average_spread, intrinsic_spread in zip(
    tenors,
    average_spreads,
    intrinsic_spreads,
):

    difference = intrinsic_spread - average_spread

    print(f"{tenor:10d}" f"{average_spread:20.6f}" f"{intrinsic_spread:20.6f}" f"{difference:20.6f}")


# ============================================================================
# 3. PLOT AVERAGE AND INTRINSIC INDEX SPREADS
# ============================================================================
#
# Plot both spread measures against maturity. This makes it easy to see
# whether the portfolio-based intrinsic spread differs materially from the
# simple average of the constituent CDS spreads.
# ============================================================================

print("\n" + "=" * 78)
print("3. AVERAGE AND INTRINSIC SPREAD TERM STRUCTURE")
print("=" * 78)

plt.figure(figsize=(9, 6))

plt.plot(
    tenors,
    average_spreads,
    marker="o",
    label="Average Spread",
)

plt.plot(
    tenors,
    intrinsic_spreads,
    marker="o",
    label="Intrinsic Spread",
)

plt.xlabel("Maturity (years)")
plt.ylabel("Spread (bp)")

plt.title("CDS Index Average and Intrinsic Spread")

plt.xticks(tenors)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 4. PLOT INTRINSIC MINUS AVERAGE SPREAD
# ============================================================================
#
# The difference isolates the effect of calculating the index spread from
# aggregate CDS portfolio valuation rather than taking a simple arithmetic
# average of the constituent spreads.
# ============================================================================

print("\n" + "=" * 78)
print("4. INTRINSIC MINUS AVERAGE SPREAD")
print("=" * 78)

spread_difference = intrinsic_spreads - average_spreads

plt.figure(figsize=(9, 6))

plt.plot(
    tenors,
    spread_difference,
    marker="o",
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.xlabel("Maturity (years)")
plt.ylabel("Intrinsic - Average Spread (bp)")

plt.title("CDS Index Spread Difference")

plt.xticks(tenors)

plt.grid(True)
plt.show()
