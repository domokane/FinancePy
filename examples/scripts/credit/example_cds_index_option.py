# ============================================================================
# FINANCEPY EXAMPLES - CDSIndexOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates the valuation of options on a CDS index using:
#
#   1. Anderson model
#   2. Adjusted Black model
#
# The constituent issuer curves are first calibrated from single-name CDS
# market spreads. They are then adjusted so that their intrinsic portfolio
# spread is consistent with the assumed CDS index spread.
#
# The example compares payer and receiver option values across:
#
#   - different index spread levels
#   - different option strikes
#
# The final plots make the dependence of the option value on the underlying
# CDS index spread visible.
# ============================================================================


import time

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.format_graphs import set_plot_style

from financepy.market.curves.cds_curve import CDSCurve

from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_index_option import CDSIndexOption
from financepy.products.credit.cds_index_portfolio import CDSIndexPortfolio


from helpers import build_ibor_curve
from helpers import load_heterogeneous_spread_curves

LINE = "=" * 100
SUBLINE = "-" * 100
set_plot_style()

# ============================================================================
# 1. CDS INDEX OPTION
# ============================================================================

print("\n" + LINE)
print("1. CDS INDEX OPTION")
print(LINE)


# ============================================================================
# 1.1 MARKET DATES
# ============================================================================

trade_dt = Date(1, 8, 2007)
step_in_dt = trade_dt.add_days(1)
value_dt = trade_dt

print(f"{'Trade Date':<30}: {trade_dt}")
print(f"{'Step-In Date':<30}: {step_in_dt}")
print(f"{'Value Date':<30}: {value_dt}")


# ============================================================================
# 1.2 BUILD INTEREST-RATE CURVE
# ============================================================================

libor_curve = build_ibor_curve(
    trade_dt,
)


# ============================================================================
# 1.3 SINGLE-NAME CDS MATURITIES
# ============================================================================

maturity_3yr = trade_dt.next_cds_date(36)
maturity_5yr = trade_dt.next_cds_date(60)
maturity_7yr = trade_dt.next_cds_date(84)
maturity_10yr = trade_dt.next_cds_date(120)


# ============================================================================
# 1.4 LOAD SINGLE-NAME CDS MARKET DATA
# ============================================================================

print("\n" + LINE)
print("2. BUILD SINGLE-NAME ISSUER CURVES")
print(LINE)

issuer_curves = load_heterogeneous_spread_curves(value_dt, step_in_dt, libor_curve)

print(f"{'Number of Issuers':<30}: {len(issuer_curves)}")


# ============================================================================
# 3. CDS INDEX MARKET DATA
# ============================================================================

print("\n" + LINE)
print("3. CDS INDEX MARKET DATA")
print(LINE)

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

print(f"{'MATURITY':>20}" f"{'UPFRONT':>20}")

print("-" * 40)

for maturity_dt, upfront in zip(
    index_maturity_dts,
    index_upfronts,
):

    print(f"{str(maturity_dt):>20}" f"{upfront:20.6f}")


# ============================================================================
# 4. CDS INDEX OPTION PARAMETERS
# ============================================================================

print("\n" + LINE)
print("4. CDS INDEX OPTION PARAMETERS")
print(LINE)

index_cpn = 0.004
volatility = 0.50

expiry_dt = Date(1, 2, 2008)
maturity_dt = Date(20, 12, 2011)

notional = 10000.0
tolerance = 1e-6

print(f"{'Index Coupon':<35}: {index_cpn * 10000.0:.6f} bp")
print(f"{'Volatility':<35}: {volatility * 100.0:.6f}%")
print(f"{'Expiry Date':<35}: {expiry_dt}")
print(f"{'Maturity Date':<35}: {maturity_dt}")
print(f"{'Notional':<35}: {notional:.2f}")
print(f"{'Recovery Rate':<35}: {index_recovery:.6f}")
print(f"{'Calibration Tolerance':<35}: {tolerance:.8f}")


# ============================================================================
# 5. VALUE CDS INDEX OPTIONS
# ============================================================================
#
# For each assumed CDS index spread:
#
#   1. Build an index CDS curve
#   2. Adjust the constituent issuer curves
#   3. Value payer and receiver options using Anderson
#   4. Value payer and receiver options using adjusted Black
#
# Payer CDS options benefit from widening credit spreads.
#
# Receiver CDS options benefit from tightening credit spreads.
# ============================================================================

print("\n" + LINE)
print("5. CDS INDEX OPTION VALUATION")
print(LINE)

index_levels_bps = [
    20,
    40,
    50,
]

strike_levels_bps = [
    20,
    60,
]

results = []

index_portfolio = CDSIndexPortfolio()


for index_bps in index_levels_bps:

    # ========================================================================
    # BUILD CDS INDEX CURVE
    # ========================================================================

    cds_contracts = []

    for index_maturity_dt in index_maturity_dts:

        cds = CDS(
            value_dt,
            index_maturity_dt,
            index_bps / 10000.0,
        )

        cds_contracts.append(cds)

    index_curve = CDSCurve(
        value_dt,
        cds_contracts,
        libor_curve,
        index_recovery,
    )

    # ========================================================================
    # ADJUST ISSUER CURVES TO THE INDEX
    # ========================================================================

    index_spreads = [index_bps / 10000.0] * len(index_maturity_dts)

    calibration_start = time.perf_counter()

    adjusted_issuer_curves = index_portfolio.hazard_rate_adjust_intrinsic(
        value_dt,
        issuer_curves,
        index_spreads,
        index_upfronts,
        index_maturity_dts,
        index_recovery,
        tolerance,
    )

    calibration_time = time.perf_counter() - calibration_start

    # ========================================================================
    # VALUE EACH STRIKE
    # ========================================================================

    for strike_bps in strike_levels_bps:

        strike = strike_bps / 10000.0

        option = CDSIndexOption(
            expiry_dt,
            maturity_dt,
            index_cpn,
            strike,
            notional,
        )

        # ====================================================================
        # ANDERSON MODEL
        # ====================================================================

        start = time.perf_counter()

        (
            anderson_payer,
            anderson_receiver,
            strike_value,
            mu,
            exp_h,
        ) = option.value_anderson(
            value_dt,
            adjusted_issuer_curves,
            index_recovery,
            volatility,
        )

        anderson_time = time.perf_counter() - start

        # ====================================================================
        # ADJUSTED BLACK MODEL
        # ====================================================================

        start = time.perf_counter()

        (
            black_payer,
            black_receiver,
        ) = option.value_adjusted_black(
            value_dt,
            index_curve,
            index_recovery,
            libor_curve,
            volatility,
        )

        black_time = time.perf_counter() - start

        # ====================================================================
        # STORE RESULTS
        # ====================================================================

        results.append(
            {
                "index_bps": index_bps,
                "strike_bps": strike_bps,
                "anderson_payer": anderson_payer,
                "anderson_receiver": anderson_receiver,
                "strike_value": strike_value,
                "mu": mu,
                "exp_h": exp_h,
                "black_payer": black_payer,
                "black_receiver": black_receiver,
                "calibration_time": calibration_time,
                "anderson_time": anderson_time,
                "black_time": black_time,
            }
        )


# ============================================================================
# 6. ANDERSON MODEL RESULTS
# ============================================================================

print("\n" + LINE)
print("6. ANDERSON MODEL RESULTS")
print(LINE)

print(
    f"{'INDEX':>10}"
    f"{'STRIKE':>10}"
    f"{'PAYER':>16}"
    f"{'RECEIVER':>16}"
    f"{'G(K)':>16}"
    f"{'MU':>16}"
    f"{'EXP(H)':>16}"
)

print("-" * 100)

for result in results:

    print(
        f"{result['index_bps']:10.2f}"
        f"{result['strike_bps']:10.2f}"
        f"{result['anderson_payer']:16.6f}"
        f"{result['anderson_receiver']:16.6f}"
        f"{result['strike_value']:16.6f}"
        f"{result['mu']:16.6f}"
        f"{result['exp_h']:16.6f}"
    )


# ============================================================================
# 7. ADJUSTED BLACK MODEL RESULTS
# ============================================================================

print("\n" + LINE)
print("7. ADJUSTED BLACK MODEL RESULTS")
print(LINE)

print(f"{'INDEX':>10}" f"{'STRIKE':>10}" f"{'PAYER':>20}" f"{'RECEIVER':>20}")

print("-" * 60)

for result in results:

    print(
        f"{result['index_bps']:10.2f}"
        f"{result['strike_bps']:10.2f}"
        f"{result['black_payer']:20.6f}"
        f"{result['black_receiver']:20.6f}"
    )


# ============================================================================
# 8. ANDERSON VERSUS ADJUSTED BLACK
# ============================================================================

print("\n" + LINE)
print("8. ANDERSON VERSUS ADJUSTED BLACK")
print(LINE)

print(
    f"{'INDEX':>8}"
    f"{'STRIKE':>8}"
    f"{'AND PAY':>14}"
    f"{'BLACK PAY':>14}"
    f"{'PAY DIFF':>14}"
    f"{'AND REC':>14}"
    f"{'BLACK REC':>14}"
    f"{'REC DIFF':>14}"
)

print("-" * 100)

for result in results:

    payer_difference = result["anderson_payer"] - result["black_payer"]

    receiver_difference = result["anderson_receiver"] - result["black_receiver"]

    print(
        f"{result['index_bps']:8.2f}"
        f"{result['strike_bps']:8.2f}"
        f"{result['anderson_payer']:14.6f}"
        f"{result['black_payer']:14.6f}"
        f"{payer_difference:14.6f}"
        f"{result['anderson_receiver']:14.6f}"
        f"{result['black_receiver']:14.6f}"
        f"{receiver_difference:14.6f}"
    )


# ============================================================================
# 9. TIMING RESULTS
# ============================================================================

print("\n" + LINE)
print("9. CALCULATION TIMES")
print(LINE)

print(f"{'INDEX':>10}" f"{'STRIKE':>10}" f"{'CALIBRATION':>18}" f"{'ANDERSON':>18}" f"{'BLACK':>18}")

print("-" * 74)

for result in results:

    print(
        f"{result['index_bps']:10.2f}"
        f"{result['strike_bps']:10.2f}"
        f"{result['calibration_time']:18.6f}"
        f"{result['anderson_time']:18.6f}"
        f"{result['black_time']:18.6f}"
    )


# ============================================================================
# 10. PLOT PAYER OPTION VALUE AGAINST INDEX SPREAD
# ============================================================================
#
# A payer CDS option benefits from spread widening. The plots below show how
# the payer option value changes as the assumed CDS index spread increases.
# ============================================================================

print("\n" + LINE)
print("10. PAYER OPTION VALUE VERSUS INDEX SPREAD")
print(LINE)

for strike_bps in strike_levels_bps:

    subset = [result for result in results if result["strike_bps"] == strike_bps]

    x = np.array([result["index_bps"] for result in subset])

    anderson_values = np.array([result["anderson_payer"] for result in subset])

    black_values = np.array([result["black_payer"] for result in subset])

    plt.figure(figsize=(9, 6))

    plt.plot(
        x,
        anderson_values,
        marker="o",
        label="Anderson",
    )

    plt.plot(
        x,
        black_values,
        marker="o",
        label="Adjusted Black",
    )

    plt.axvline(
        strike_bps,
        linestyle="--",
        label="Strike",
    )

    plt.xlabel("CDS Index Spread (bp)")
    plt.ylabel("Payer Option Value")

    plt.title(f"CDS Index Payer Option - Strike {strike_bps} bp")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================
# 11. PLOT RECEIVER OPTION VALUE AGAINST INDEX SPREAD
# ============================================================================
#
# A receiver CDS option benefits from spread tightening. Its value therefore
# generally behaves in the opposite direction to the payer option as the
# underlying index spread changes.
# ============================================================================

print("\n" + LINE)
print("11. RECEIVER OPTION VALUE VERSUS INDEX SPREAD")
print(LINE)

for strike_bps in strike_levels_bps:

    subset = [result for result in results if result["strike_bps"] == strike_bps]

    x = np.array([result["index_bps"] for result in subset])

    anderson_values = np.array([result["anderson_receiver"] for result in subset])

    black_values = np.array([result["black_receiver"] for result in subset])

    plt.figure(figsize=(9, 6))

    plt.plot(
        x,
        anderson_values,
        marker="o",
        label="Anderson",
    )

    plt.plot(
        x,
        black_values,
        marker="o",
        label="Adjusted Black",
    )

    plt.axvline(
        strike_bps,
        linestyle="--",
        label="Strike",
    )

    plt.xlabel("CDS Index Spread (bp)")
    plt.ylabel("Receiver Option Value")

    plt.title(f"CDS Index Receiver Option - Strike {strike_bps} bp")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================
# 12. PLOT MODEL DIFFERENCE
# ============================================================================
#
# This plot isolates the difference between the Anderson and adjusted Black
# valuations:
#
#       model difference = Anderson value - Adjusted Black value
#
# A value near zero indicates close agreement between the two approaches.
# ============================================================================

print("\n" + LINE)
print("12. ANDERSON MINUS ADJUSTED BLACK")
print(LINE)

for strike_bps in strike_levels_bps:

    subset = [result for result in results if result["strike_bps"] == strike_bps]

    x = np.array([result["index_bps"] for result in subset])

    payer_difference = np.array([result["anderson_payer"] - result["black_payer"] for result in subset])

    receiver_difference = np.array([result["anderson_receiver"] - result["black_receiver"] for result in subset])

    plt.figure(figsize=(9, 6))

    plt.plot(
        x,
        payer_difference,
        marker="o",
        label="Payer Difference",
    )

    plt.plot(
        x,
        receiver_difference,
        marker="o",
        label="Receiver Difference",
    )

    plt.axhline(
        0.0,
        linestyle="--",
    )

    plt.xlabel("CDS Index Spread (bp)")
    plt.ylabel("Anderson - Adjusted Black")

    plt.title(f"CDS Index Option Model Difference - Strike {strike_bps} bp")

    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================
# 13. SUMMARY
# ============================================================================

print("\n" + LINE)
print("13. SUMMARY")
print(LINE)

print("The CDS index option is valued using both the Anderson model " "and the adjusted Black model.")

print(
    "Before the Anderson valuation, the constituent issuer curves are "
    "adjusted so that their intrinsic portfolio spread is consistent "
    "with the assumed CDS index spread."
)

print(
    "The payer option provides exposure to widening CDS index spreads, "
    "while the receiver option provides exposure to tightening spreads."
)

print(
    "The plots show how payer and receiver option values change with "
    "the underlying CDS index spread and how the two pricing methods "
    "compare."
)

print("\n" + LINE)
print("END OF CDS INDEX OPTION EXAMPLE")
print(LINE)
