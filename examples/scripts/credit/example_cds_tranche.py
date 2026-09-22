# ============================================================================
# FINANCEPY EXAMPLES - CDSTranche
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of synthetic CDS index tranches.
#
# A CDS tranche absorbs losses only over a specified interval of cumulative
# portfolio loss:
#
#       [attachment point, detachment point]
#
# For example:
#
#       0%-3%     Equity tranche
#       3%-6%     Junior mezzanine tranche
#       6%-9%     Mezzanine tranche
#       9%-12%    Senior mezzanine tranche
#       12%-22%   Senior tranche
#       22%-60%   Super-senior tranche
#
# The 0%-3% tranche begins absorbing losses immediately when defaults occur.
# Once cumulative portfolio losses reach 3%, the tranche is exhausted.
#
# The 3%-6% tranche is protected from the first 3% of portfolio losses but
# begins taking losses once cumulative portfolio loss exceeds 3%.
#
# Moving higher through the capital structure therefore increases the amount
# of portfolio loss that must occur before the tranche begins to suffer loss.
#
# This example considers:
#
#   1. Homogeneous issuer credit curves.
#   2. Heterogeneous issuer credit curves.
#   3. Alternative loss-distribution construction methods.
#   4. Different attachment and detachment points.
#   5. The effect of default correlation on tranche spreads.
#
# Correlation is particularly important for tranche valuation because it
# changes the distribution of portfolio losses. It can therefore affect
# different parts of the capital structure in very different ways.
#
# Date convention:
#
#       value_dt   = trade_dt
#       step_in_dt = trade_dt + 1 day
#
# ============================================================================

import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.products.credit.cds_tranche import (
    CDSTranche,
    FinLossDistributionBuilder,
)
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio
from financepy.utils.format_graphs import set_plot_style

from helpers import build_ibor_curve
from helpers import build_homogeneous_issuer_curves
from helpers import load_heterogeneous_issuer_curves

# ============================================================================
# GLOBAL OUTPUT SETTINGS
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100
set_plot_style()


# ============================================================================
# 1. CDS TRANCHE
# ============================================================================
#
# Construct a collection of tranches covering different regions of the
# portfolio-loss distribution.
#
# A tranche with attachment k1 and detachment k2 is exposed to portfolio
# losses between k1 and k2.
#
# The final 0%-60% tranche spans the complete loss interval covered by the
# other tranches. It should therefore not be interpreted as being senior to
# the 22%-60% tranche. It is an overlapping, index-like tranche useful for
# comparison.
# ============================================================================

print("\n" + LINE)
print("1. CDS TRANCHE")
print(LINE)


# ============================================================================
# 1.1 DATES
# ============================================================================

trade_dt = Date(
    1,
    3,
    2007,
)

value_dt = trade_dt

step_in_dt = trade_dt.add_days(
    1,
)

tranche_maturity = Date(
    20,
    12,
    2011,
)

print(f"{'Trade Date':<35}: {trade_dt}")
print(f"{'Value Date':<35}: {value_dt}")
print(f"{'Step-In Date':<35}: {step_in_dt}")
print(f"{'Tranche Maturity':<35}: {tranche_maturity}")


# ============================================================================
# 1.2 BUILD INTEREST-RATE CURVE
# ============================================================================

libor_curve = build_ibor_curve(
    value_dt,
)


# ============================================================================
# 1.3 DEFINE TRANCHES
# ============================================================================

tranche1 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.00,
    0.03,
)

tranche2 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.03,
    0.06,
)

tranche3 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.06,
    0.09,
)

tranche4 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.09,
    0.12,
)

tranche5 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.12,
    0.22,
)

tranche6 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.22,
    0.60,
)

tranche7 = CDSTranche(
    value_dt,
    tranche_maturity,
    0.00,
    0.60,
)

tranches = [
    tranche1,
    tranche2,
    tranche3,
    tranche4,
    tranche5,
    tranche6,
    tranche7,
]


# ============================================================================
# 1.4 TRANCHE LABELS
# ============================================================================

tranche_labels = [f"{100.0 * tranche.k1:.0f}-{100.0 * tranche.k2:.0f}%" for tranche in tranches]


# ============================================================================
# 1.5 MODEL PARAMETERS
# ============================================================================
#
# corr1 and corr2 are the correlation parameters passed to the tranche
# valuation model.
#
# The initial example uses different values for the lower and upper tranche
# boundaries.
#
# Later we set corr1 = corr2 and vary the common correlation to examine
# correlation sensitivity.
# ============================================================================

corr1 = 0.30
corr2 = 0.35

upfront = 0.0
spd = 0.0

num_points = 40

cds_index = CDSIndexPortfolio()

print("\n" + SUBLINE)

print(f"{'Correlation 1':<35}: {corr1:12.6f}")
print(f"{'Correlation 2':<35}: {corr2:12.6f}")
print(f"{'Upfront':<35}: {upfront:12.6f}")
print(f"{'Running Spread Input':<35}: {spd:12.6f}")
print(f"{'Loss Distribution Points':<35}: {num_points:d}")


# ============================================================================
# 2. HOMOGENEOUS CREDIT PORTFOLIO
# ============================================================================
#
# A homogeneous portfolio assumes every issuer has the same CDS spread curve.
#
# This is useful because it removes cross-sectional differences between
# issuers. Differences between tranche spreads then arise primarily from:
#
#   - attachment and detachment points,
#   - portfolio loss mechanics,
#   - correlation,
#   - and the loss-distribution approximation.
#
# The example uses 125 identical credits.
# ============================================================================

print("\n" + LINE)
print("2. HOMOGENEOUS CREDIT PORTFOLIO")
print(LINE)

num_credits = 125

spd_3yr = 0.0012
spd_5yr = 0.0025
spd_7yr = 0.0034
spd_10yr = 0.0046

issuer_curves = build_homogeneous_issuer_curves(
    value_dt,
    step_in_dt,
    libor_curve,
    spd_3yr,
    spd_5yr,
    spd_7yr,
    spd_10yr,
    num_credits,
)

print(f"{'Number of Credits':<35}: {len(issuer_curves)}")

print(f"{'3Y CDS Spread':<35}: {spd_3yr * 10000.0:12.6f} bp")
print(f"{'5Y CDS Spread':<35}: {spd_5yr * 10000.0:12.6f} bp")
print(f"{'7Y CDS Spread':<35}: {spd_7yr * 10000.0:12.6f} bp")
print(f"{'10Y CDS Spread':<35}: {spd_10yr * 10000.0:12.6f} bp")


# ============================================================================
# 2.1 INTRINSIC INDEX SPREAD
# ============================================================================
#
# The intrinsic spread is obtained by aggregating the protection and premium
# legs of the constituent CDS portfolio.
#
# This gives a useful reference spread for the underlying index portfolio
# before splitting the portfolio into individual loss tranches.
# ============================================================================

intrinsic_spd = (
    cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        tranche_maturity,
        issuer_curves,
    )
    * 10000.0
)

adjusted_spd = intrinsic_spd / 0.60

print("\n" + SUBLINE)

print(f"{'Intrinsic Spread at Tranche Maturity':<45}: " f"{intrinsic_spd:14.6f} bp")

print(f"{'Intrinsic Spread / 60%':<45}: " f"{adjusted_spd:14.6f} bp")


# ============================================================================
# 2.2 VALUE HOMOGENEOUS TRANCHES
# ============================================================================
#
# value_bc() returns the tranche valuation results. The fourth element of the
# result, v[3], is the model tranche spread.
#
# Multiplying by 10,000 converts the decimal spread into basis points.
#
# Results are stored so they can subsequently be plotted.
# ============================================================================

print("\n" + SUBLINE)
print("HOMOGENEOUS TRANCHE RESULTS")
print(SUBLINE)

print(f"{'METHOD':<35}" f"{'TIME':>12}" f"{'POINTS':>10}" f"{'K1':>10}" f"{'K2':>10}" f"{'SPREAD (bp)':>18}")

print("-" * 95)

homogeneous_results = {}

for method in FinLossDistributionBuilder:

    method_results = []

    for tranche in tranches:

        start = time.time()

        value = tranche.value_bc(
            value_dt,
            issuer_curves,
            upfront,
            spd,
            corr1,
            corr2,
            num_points,
            method,
        )

        elapsed = time.time() - start

        tranche_spread = value[3] * 10000.0

        method_results.append(
            tranche_spread,
        )

        print(
            f"{str(method):<35}"
            f"{elapsed:12.6f}"
            f"{num_points:10d}"
            f"{tranche.k1:10.4f}"
            f"{tranche.k2:10.4f}"
            f"{tranche_spread:18.6f}"
        )

    homogeneous_results[method] = np.asarray(
        method_results,
    )


# ============================================================================
# 3. HETEROGENEOUS CREDIT PORTFOLIO
# ============================================================================
#
# The heterogeneous portfolio uses the individual issuer CDS curves contained
# in the market-data set.
#
# Unlike the homogeneous example, constituents can therefore have different
# credit spreads and recovery assumptions.
#
# This introduces cross-sectional credit dispersion into the portfolio loss
# distribution.
# ============================================================================

print("\n" + LINE)
print("3. HETEROGENEOUS CREDIT PORTFOLIO")
print(LINE)

heterogeneous_issuer_curves = load_heterogeneous_issuer_curves(
    value_dt,
    step_in_dt,
    libor_curve,
)

print(f"{'Number of Credits':<35}: " f"{len(heterogeneous_issuer_curves)}")


# ============================================================================
# 3.1 INTRINSIC INDEX SPREAD
# ============================================================================

heterogeneous_intrinsic_spd = (
    cds_index.intrinsic_spread(
        value_dt,
        step_in_dt,
        tranche_maturity,
        heterogeneous_issuer_curves,
    )
    * 10000.0
)

heterogeneous_adjusted_spd = heterogeneous_intrinsic_spd / 0.60

print("\n" + SUBLINE)

print(f"{'Intrinsic Spread at Tranche Maturity':<45}: " f"{heterogeneous_intrinsic_spd:14.6f} bp")

print(f"{'Intrinsic Spread / 60%':<45}: " f"{heterogeneous_adjusted_spd:14.6f} bp")


# ============================================================================
# 3.2 VALUE HETEROGENEOUS TRANCHES
# ============================================================================

print("\n" + SUBLINE)
print("HETEROGENEOUS TRANCHE RESULTS")
print(SUBLINE)

print(f"{'METHOD':<35}" f"{'TIME':>12}" f"{'POINTS':>10}" f"{'K1':>10}" f"{'K2':>10}" f"{'SPREAD (bp)':>18}")

print("-" * 95)

heterogeneous_results = {}

for method in FinLossDistributionBuilder:

    method_results = []

    for tranche in tranches:

        start = time.time()

        value = tranche.value_bc(
            value_dt,
            heterogeneous_issuer_curves,
            upfront,
            spd,
            corr1,
            corr2,
            num_points,
            method,
        )

        elapsed = time.time() - start

        tranche_spread = value[3] * 10000.0

        method_results.append(
            tranche_spread,
        )

        print(
            f"{str(method):<35}"
            f"{elapsed:12.6f}"
            f"{num_points:10d}"
            f"{tranche.k1:10.4f}"
            f"{tranche.k2:10.4f}"
            f"{tranche_spread:18.6f}"
        )

    heterogeneous_results[method] = np.asarray(
        method_results,
    )


# ============================================================================
# 4. HOMOGENEOUS TRANCHE SPREADS
# ============================================================================
#
# Plot tranche spreads across the capital structure for each loss-distribution
# method.
#
# Moving from the 0%-3% tranche toward the more senior tranches increases the
# amount of portfolio loss that must occur before the tranche is affected.
#
# The 0%-60% tranche is included as an index-like comparison. Because it
# overlaps the other tranches, it should not be interpreted as the final
# sequential tranche in the capital structure.
# ============================================================================

print("\n" + LINE)
print("4. HOMOGENEOUS TRANCHE SPREADS")
print(LINE)

for method, spreads in homogeneous_results.items():

    plt.figure(
        figsize=(10, 6),
    )

    plt.plot(
        tranche_labels,
        spreads,
        marker="o",
    )

    plt.xlabel(
        "Tranche Attachment-Detachment Interval",
    )

    plt.ylabel(
        "Tranche Spread (bp)",
    )

    plt.title("Homogeneous Portfolio Tranche Spreads\n" + str(method))

    plt.grid(
        True,
    )

    plt.tight_layout()

    plt.show()


# ============================================================================
# 5. HETEROGENEOUS TRANCHE SPREADS
# ============================================================================
#
# Repeat the capital-structure plot using the heterogeneous constituent
# credit curves.
#
# Differences from the homogeneous case arise because individual issuers no
# longer share identical credit curves.
# ============================================================================

print("\n" + LINE)
print("5. HETEROGENEOUS TRANCHE SPREADS")
print(LINE)

for method, spreads in heterogeneous_results.items():

    plt.figure(
        figsize=(10, 6),
    )

    plt.plot(
        tranche_labels,
        spreads,
        marker="o",
    )

    plt.xlabel(
        "Tranche Attachment-Detachment Interval",
    )

    plt.ylabel(
        "Tranche Spread (bp)",
    )

    plt.title("Heterogeneous Portfolio Tranche Spreads\n" + str(method))

    plt.grid(
        True,
    )

    plt.tight_layout()

    plt.show()


# ============================================================================
# 6. HOMOGENEOUS VERSUS HETEROGENEOUS PORTFOLIOS
# ============================================================================
#
# Plot the two portfolio assumptions on the same graph.
#
# This isolates the effect of constituent heterogeneity.
#
# If the two curves are close for a particular tranche, replacing the
# individual issuer curves with a homogeneous approximation has relatively
# little effect on the calculated spread for that tranche under the chosen
# assumptions.
#
# A larger separation indicates greater sensitivity to the cross-sectional
# distribution of issuer credit quality.
# ============================================================================

print("\n" + LINE)
print("6. HOMOGENEOUS VERSUS HETEROGENEOUS PORTFOLIOS")
print(LINE)

for method in FinLossDistributionBuilder:

    plt.figure(
        figsize=(10, 6),
    )

    plt.plot(
        tranche_labels,
        homogeneous_results[method],
        marker="o",
        label="Homogeneous",
    )

    plt.plot(
        tranche_labels,
        heterogeneous_results[method],
        marker="o",
        label="Heterogeneous",
    )

    plt.xlabel(
        "Tranche Attachment-Detachment Interval",
    )

    plt.ylabel(
        "Tranche Spread (bp)",
    )

    plt.title("Homogeneous versus Heterogeneous Portfolio\n" + str(method))

    plt.grid(
        True,
    )

    plt.legend()

    plt.tight_layout()

    plt.show()


# ============================================================================
# 7. HETEROGENEITY EFFECT
# ============================================================================
#
# Plot:
#
#       heterogeneous spread - homogeneous spread
#
# rather than the absolute spreads.
#
# This makes the effect of constituent heterogeneity easier to identify when
# the absolute tranche spreads are very different in magnitude.
#
# A positive value means that the heterogeneous portfolio produces a larger
# tranche spread under the chosen model assumptions.
#
# A negative value means that the homogeneous portfolio produces the larger
# spread.
# ============================================================================

print("\n" + LINE)
print("7. HETEROGENEITY EFFECT")
print(LINE)

for method in FinLossDistributionBuilder:

    spread_difference = heterogeneous_results[method] - homogeneous_results[method]

    print("\n" + str(method))

    print(f"{'TRANCHE':>15}" f"{'HOMOGENEOUS':>20}" f"{'HETEROGENEOUS':>20}" f"{'DIFFERENCE':>20}")

    print("-" * 75)

    for (
        label,
        homogeneous_spread,
        heterogeneous_spread,
        difference,
    ) in zip(
        tranche_labels,
        homogeneous_results[method],
        heterogeneous_results[method],
        spread_difference,
    ):

        print(f"{label:>15}" f"{homogeneous_spread:20.6f}" f"{heterogeneous_spread:20.6f}" f"{difference:20.6f}")

    plt.figure(
        figsize=(10, 6),
    )

    plt.plot(
        tranche_labels,
        spread_difference,
        marker="o",
    )

    plt.axhline(
        0.0,
        linestyle="--",
    )

    plt.xlabel(
        "Tranche Attachment-Detachment Interval",
    )

    plt.ylabel(
        "Heterogeneous - Homogeneous Spread (bp)",
    )

    plt.title("Effect of Issuer Heterogeneity\n" + str(method))

    plt.grid(
        True,
    )

    plt.tight_layout()

    plt.show()


# ============================================================================
# 8. CORRELATION SENSITIVITY
# ============================================================================
#
# Default correlation is one of the central risk factors in tranche pricing.
#
# Correlation changes the distribution of portfolio losses.
#
# With relatively low correlation, issuer defaults are more independent and
# portfolio losses tend to be more diversified.
#
# With higher correlation, common systematic outcomes become more important.
# The probability distribution places relatively more weight on scenarios in
# which many issuers experience similar credit outcomes.
#
# Importantly, this does NOT imply that increasing correlation must increase
# every tranche spread.
#
# Equity, mezzanine and senior tranches occupy different regions of the loss
# distribution. A change in correlation can therefore affect them in
# different directions.
#
# The following experiment sets:
#
#       corr1 = corr2 = correlation
#
# and varies that common correlation while keeping the issuer credit curves
# unchanged.
#
# The heterogeneous portfolio is used for this experiment.
# ============================================================================

print("\n" + LINE)
print("8. CORRELATION SENSITIVITY")
print(LINE)

correlations = np.linspace(
    0.05,
    0.80,
    16,
)

methods = list(
    FinLossDistributionBuilder,
)

correlation_method = methods[0]

print(f"{'Loss Distribution Method':<35}: " f"{correlation_method}")

correlation_results = {label: [] for label in tranche_labels}

print("\n" + SUBLINE)

print(f"{'CORRELATION':>15}" f"{'TRANCHE':>15}" f"{'SPREAD (bp)':>20}")

print("-" * 50)

for correlation in correlations:

    for label, tranche in zip(
        tranche_labels,
        tranches,
    ):

        value = tranche.value_bc(
            value_dt,
            heterogeneous_issuer_curves,
            upfront,
            spd,
            correlation,
            correlation,
            num_points,
            correlation_method,
        )

        tranche_spread = value[3] * 10000.0

        correlation_results[label].append(
            tranche_spread,
        )

        print(f"{correlation:15.4f}" f"{label:>15}" f"{tranche_spread:20.6f}")


# ============================================================================
# 9. PLOT CORRELATION SENSITIVITY
# ============================================================================
#
# Each line represents one tranche.
#
# The graph demonstrates that correlation is not simply a parallel spread
# shift across the capital structure. Each attachment-detachment interval
# responds according to its exposure to the changing portfolio loss
# distribution.
# ============================================================================

print("\n" + LINE)
print("9. TRANCHE SPREAD VERSUS DEFAULT CORRELATION")
print(LINE)

plt.figure(
    figsize=(11, 7),
)

for label in tranche_labels:

    plt.plot(
        correlations,
        correlation_results[label],
        marker="o",
        label=label,
    )

plt.xlabel(
    "Default Correlation",
)

plt.ylabel(
    "Tranche Spread (bp)",
)

plt.title("CDS Tranche Spread Sensitivity to Default Correlation")

plt.grid(
    True,
)

plt.legend(
    title="Tranche",
)

plt.tight_layout()

plt.show()


# ============================================================================
# 10. SELECTED TRANCHE CORRELATION SENSITIVITY
# ============================================================================
#
# The previous graph contains all tranches and can become visually dominated
# by the equity tranche because tranche spreads can differ greatly in scale.
#
# Plotting selected tranches separately makes the correlation behaviour of
# the mezzanine and senior portions of the capital structure easier to see.
#
# The 0%-60% overlapping tranche is omitted here because it represents a broad
# portfolio-loss interval rather than one sequential layer of the capital
# structure.
# ============================================================================

print("\n" + LINE)
print("10. SELECTED TRANCHE CORRELATION SENSITIVITY")
print(LINE)

selected_labels = [
    "3-6%",
    "6-9%",
    "9-12%",
    "12-22%",
    "22-60%",
]

plt.figure(
    figsize=(11, 7),
)

for label in selected_labels:

    if label in correlation_results:

        plt.plot(
            correlations,
            correlation_results[label],
            marker="o",
            label=label,
        )

plt.xlabel(
    "Default Correlation",
)

plt.ylabel(
    "Tranche Spread (bp)",
)

plt.title("Mezzanine and Senior Tranche Correlation Sensitivity")

plt.grid(
    True,
)

plt.legend(
    title="Tranche",
)

plt.tight_layout()

plt.show()


# ============================================================================
# 11. SUMMARY
# ============================================================================
#
# The example illustrates three important features of CDS tranche valuation:
#
# 1. SUBORDINATION
#
#    Attachment and detachment points determine which region of portfolio
#    losses is allocated to each tranche. Junior tranches absorb losses before
#    more senior tranches.
#
# 2. ISSUER HETEROGENEITY
#
#    Replacing individual issuer credit curves with one homogeneous curve can
#    change the portfolio loss distribution and therefore tranche spreads.
#    The homogeneous-versus-heterogeneous plots show where this approximation
#    matters most under the chosen market inputs.
#
# 3. DEFAULT CORRELATION
#
#    Correlation changes the shape of the portfolio loss distribution.
#    Because different tranches occupy different regions of that distribution,
#    their sensitivities to correlation can be substantially different.
#
# The resulting tranche spread is therefore not determined solely by the
# average CDS spread of the underlying portfolio. It also depends on the
# distribution of constituent credit risk, the tranche attachment and
# detachment points, recovery assumptions, correlation assumptions and the
# loss-distribution model.
# ============================================================================

print("\n" + LINE)
print("11. SUMMARY")
print(LINE)

print(
    "The tranche calculations demonstrate how portfolio credit risk is "
    "redistributed across different attachment and detachment intervals."
)

print("Compare the homogeneous and heterogeneous plots to see the effect " "of constituent credit dispersion.")

print(
    "Compare the correlation-sensitivity plots to see how changes in the "
    "portfolio loss distribution affect different parts of the capital "
    "structure."
)

print(LINE)
