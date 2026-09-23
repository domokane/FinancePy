# ============================================================================
# FINANCEPY EXAMPLES - CDSIndex
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.global_vars import CLEAN, DIRTY
from financepy.utils.date import Date
from financepy.utils.global_vars import ONE_MILLION
from financepy.utils.format_graphs import set_plot_style

from financepy.products.credit.cds import CDS

from helpers import build_ibor_curve
from helpers import build_issuer_curve


# ============================================================================
# GRAPH FORMATTING
# ============================================================================

set_plot_style()

LINE = "=" * 78


# ============================================================================
# 1. VALUE CDS INDEX
# ============================================================================
# What this section demonstrates:
#
# Values a CDS index using a discount curve and an issuer credit curve.
#
# The CDS valuation is decomposed into:
#
#   - par spread
#   - dirty value
#   - clean value
#   - clean price
#   - accrued premium
#   - protection-leg PV
#   - premium-leg PV
#   - risky PV01
#
# The example uses:
#
#     value date = trade date
#
# while the CDS protection starts on the step-in date.
# ============================================================================

print("\n" + LINE)
print("1. VALUE CDS INDEX")
print(LINE)

trade_dt = Date(7, 2, 2006)
value_dt = trade_dt
step_in_dt = trade_dt.add_days(1)

maturity_dt = Date(20, 6, 2010)

libor_curve = build_ibor_curve(
    value_dt,
)

cds_recovery = 0.40

issuer_curve = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
)

notional = 10.0 * ONE_MILLION
long_protection = True
index_cpn = 0.004

cds_index_contract = CDS(
    step_in_dt,
    maturity_dt,
    index_cpn,
    notional,
    long_protection,
)

print(f"{'Trade Date':<30}: {trade_dt}")
print(f"{'Value Date':<30}: {value_dt}")
print(f"{'Step-In Date':<30}: {step_in_dt}")
print(f"{'Maturity Date':<30}: {maturity_dt}")
print(f"{'Notional':<30}: {notional:,.2f}")
print(f"{'Index Coupon':<30}: {index_cpn * 10000.0:.2f} bp")
print(f"{'Recovery Rate':<30}: {cds_recovery:.2%}")

par_spread = (
    cds_index_contract.par_spread(
        value_dt,
        issuer_curve,
        cds_recovery,
    )
    * 10000.0
)

values = cds_index_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery,
)

dirty_value = values[DIRTY]
clean_value = values[CLEAN]

clean_price = cds_index_contract.clean_price(
    value_dt,
    issuer_curve,
    cds_recovery,
)

accrued_days = cds_index_contract.accrued_days(
    value_dt,
)

accrued_interest = cds_index_contract.accrued_interest(
    value_dt,
)

protection_pv = cds_index_contract.prot_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

premium_pv = cds_index_contract.premium_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

dirty_rpv01, clean_rpv01 = cds_index_contract.rpv01(
    value_dt,
    issuer_curve,
)

print("\n" + "-" * 60)

print(f"{'Par Spread':<30}: {par_spread:15.6f} bp")
print(f"{'Dirty Value':<30}: {dirty_value:15.6f}")
print(f"{'Clean Value':<30}: {clean_value:15.6f}")
print(f"{'Clean Price':<30}: {clean_price:15.6f}")
print(f"{'Accrued Days':<30}: {accrued_days}")
print(f"{'Accrued Coupon':<30}: {accrued_interest:15.6f}")
print(f"{'Protection Leg PV':<30}: {protection_pv:15.6f}")
print(f"{'Premium Leg PV':<30}: {premium_pv:15.6f}")
print(f"{'Dirty RPV01':<30}: {dirty_rpv01:15.8f}")
print(f"{'Clean RPV01':<30}: {clean_rpv01:15.8f}")


# ============================================================================
# 2. CLEAN AND DIRTY VALUE
# ============================================================================
# What this section demonstrates:
#
# The dirty value contains accrued premium whereas the clean value removes
# this accrued amount.
#
# Comparing the two values provides a useful check on the treatment of
# accrued coupon in the CDS valuation.
# ============================================================================

print("\n" + LINE)
print("2. CLEAN AND DIRTY VALUE")
print(LINE)

clean_dirty_difference = (
    dirty_value
    - clean_value
)

print(
    f"{'Dirty Value':<35}"
    f"{dirty_value:20.6f}"
)

print(
    f"{'Clean Value':<35}"
    f"{clean_value:20.6f}"
)

print(
    f"{'Dirty - Clean':<35}"
    f"{clean_dirty_difference:20.6f}"
)

print(
    f"{'Accrued Coupon':<35}"
    f"{accrued_interest:20.6f}"
)


# ============================================================================
# 3. PROTECTION AND PREMIUM LEGS
# ============================================================================
# What this section demonstrates:
#
# A CDS is composed of two legs.
#
# The protection leg represents the expected payment following default.
#
# The premium leg represents the expected stream of contractual coupon
# payments.
#
# Their relative magnitude determines the value of the CDS contract.
# ============================================================================

print("\n" + LINE)
print("3. PROTECTION AND PREMIUM LEGS")
print(LINE)

print(
    f"{'Protection Leg PV':<35}"
    f"{protection_pv:20.6f}"
)

print(
    f"{'Premium Leg PV':<35}"
    f"{premium_pv:20.6f}"
)

print(
    f"{'Protection - Premium':<35}"
    f"{protection_pv - premium_pv:20.6f}"
)


# ============================================================================
# 4. INDEX COUPON VERSUS PAR SPREAD
# ============================================================================
# What this section demonstrates:
#
# The contractual coupon is fixed when the CDS index contract is created.
#
# The par spread is the spread that would make a newly entered CDS have zero
# value under the current discount and credit curves.
#
# If the contractual coupon differs from the par spread, the existing CDS
# generally has a non-zero mark-to-market value.
# ============================================================================

print("\n" + LINE)
print("4. INDEX COUPON VERSUS PAR SPREAD")
print(LINE)

print(
    f"{'Index Coupon':<35}"
    f"{index_cpn * 10000.0:20.6f} bp"
)

print(
    f"{'Par Spread':<35}"
    f"{par_spread:20.6f} bp"
)

print(
    f"{'Par Spread - Coupon':<35}"
    f"{par_spread - index_cpn * 10000.0:20.6f} bp"
)


# ============================================================================
# 5. CDS INDEX VALUE VERSUS CONTRACT COUPON
# ============================================================================
# What this section demonstrates:
#
# Holds the market discount and issuer credit curves fixed while changing
# only the contractual coupon of the CDS.
#
# For a long-protection position, increasing the contractual coupon makes the
# premium leg more expensive and therefore reduces the value of the contract.
#
# The value should be close to zero when the contractual coupon is close to
# the current par spread.
# ============================================================================

print("\n" + LINE)
print("5. CDS INDEX VALUE VERSUS CONTRACT COUPON")
print(LINE)

coupon_grid_bp = np.linspace(
    10.0,
    100.0,
    46,
)

coupon_values = []

for coupon_bp in coupon_grid_bp:

    coupon = coupon_bp / 10000.0

    contract = CDS(
        step_in_dt,
        maturity_dt,
        coupon,
        notional,
        long_protection,
    )

    contract_value = contract.value(
        value_dt,
        issuer_curve,
        cds_recovery,
    )[CLEAN]

    coupon_values.append(
        contract_value,
    )

coupon_values = np.array(
    coupon_values,
)

print(
    f"{'COUPON (bp)':>15}"
    f"{'CLEAN VALUE':>20}"
)

print("-" * 35)

for coupon_bp, contract_value in zip(
    coupon_grid_bp,
    coupon_values,
):

    print(
        f"{coupon_bp:15.2f}"
        f"{contract_value:20.6f}"
    )


# ============================================================================
# 6. PLOT CDS INDEX VALUE VERSUS CONTRACT COUPON
# ============================================================================

print("\n" + LINE)
print("6. CDS INDEX VALUE VERSUS CONTRACT COUPON")
print(LINE)

plt.figure()

plt.plot(
    coupon_grid_bp,
    coupon_values,
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.axvline(
    par_spread,
    linestyle="--",
    label="Par Spread",
)

plt.xlabel("Contract Coupon (bp)")
plt.ylabel("Clean Value")

plt.title(
    "CDS Index Value versus Contract Coupon"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 7. VALUE AT THE PAR SPREAD
# ============================================================================
# What this section demonstrates:
#
# A newly entered CDS with a contractual coupon equal to its par spread
# should have approximately zero clean value.
#
# This is an important internal consistency check because the par spread is
# defined as the spread that balances the protection and premium legs.
# ============================================================================

print("\n" + LINE)
print("7. VALUE AT THE PAR SPREAD")
print(LINE)

par_contract = CDS(
    step_in_dt,
    maturity_dt,
    par_spread / 10000.0,
    notional,
    long_protection,
)

par_values = par_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery,
)

par_dirty_value = par_values[DIRTY]
par_clean_value = par_values[CLEAN]

print(
    f"{'Par Coupon':<35}"
    f"{par_spread:20.8f} bp"
)

print(
    f"{'Dirty Value':<35}"
    f"{par_dirty_value:20.8f}"
)

print(
    f"{'Clean Value':<35}"
    f"{par_clean_value:20.8f}"
)


# ============================================================================
# 8. PAR-SPREAD CONSISTENCY CHECK
# ============================================================================
# What this section demonstrates:
#
# Tests numerically that a CDS struck at its calculated par spread has a
# clean value close to zero.
# ============================================================================

print("\n" + LINE)
print("8. PAR-SPREAD CONSISTENCY CHECK")
print(LINE)

par_value_tolerance = 1.0e-6 * notional

print(
    f"{'Absolute Clean Value':<35}"
    f"{abs(par_clean_value):20.10f}"
)

print(
    f"{'Tolerance':<35}"
    f"{par_value_tolerance:20.10f}"
)

if abs(par_clean_value) <= par_value_tolerance:

    print("CHECK: PASSED")

else:

    print("CHECK: FAILED")


# ============================================================================
# 9. CDS INDEX VALUE VERSUS ISSUER SPREAD
# ============================================================================
# What this section demonstrates:
#
# Rebuilds the issuer credit curve after applying a parallel bump to its
# calibration spreads.
#
# This is different from changing the contractual coupon in Section 5.
#
# Here the CDS contract itself is unchanged. Instead, the market credit curve
# used to value the contract is changed.
#
# For a long-protection position, wider issuer spreads generally increase the
# value of protection because the implied default risk is higher.
# ============================================================================

print("\n" + LINE)
print("9. CDS INDEX VALUE VERSUS ISSUER SPREAD")
print(LINE)

spread_bumps_bp = np.arange(
    -20.0,
    21.0,
    2.0,
)

spread_values = []
valid_spread_bumps_bp = []

for bump_bp in spread_bumps_bp:

    spread_bump = bump_bp / 10000.0

    try:

        bumped_issuer_curve = build_issuer_curve(
            value_dt,
            step_in_dt,
            libor_curve,
            cds_recovery,
            spread_bump=spread_bump,
        )

        bumped_value = cds_index_contract.value(
            value_dt,
            bumped_issuer_curve,
            cds_recovery,
        )[CLEAN]

        valid_spread_bumps_bp.append(
            bump_bp,
        )

        spread_values.append(
            bumped_value,
        )

    except Exception:

        # Large negative spread bumps can occasionally make a calibrated
        # survival curve invalid. Such points are omitted from the graph.
        pass

valid_spread_bumps_bp = np.array(
    valid_spread_bumps_bp,
)

spread_values = np.array(
    spread_values,
)

print(
    f"{'SPREAD BUMP (bp)':>20}"
    f"{'CLEAN VALUE':>20}"
)

print("-" * 40)

for bump_bp, bumped_value in zip(
    valid_spread_bumps_bp,
    spread_values,
):

    print(
        f"{bump_bp:20.2f}"
        f"{bumped_value:20.6f}"
    )


# ============================================================================
# 10. PLOT CDS INDEX VALUE VERSUS ISSUER SPREAD
# ============================================================================

print("\n" + LINE)
print("10. CDS INDEX VALUE VERSUS ISSUER SPREAD")
print(LINE)

plt.figure()

plt.plot(
    valid_spread_bumps_bp,
    spread_values,
    marker="o",
)

plt.axvline(
    0.0,
    linestyle="--",
)

plt.axhline(
    clean_value,
    linestyle="--",
    label="Base Value",
)

plt.xlabel("Parallel Issuer Spread Bump (bp)")
plt.ylabel("Clean Value")

plt.title(
    "CDS Index Value versus Issuer Spread"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 11. SPREAD 01 BY BUMP AND REVALUE
# ============================================================================
# What this section demonstrates:
#
# Spread 01 measures the change in CDS value resulting from a small change in
# the issuer credit spreads used to calibrate the credit curve.
#
# A central difference is used:
#
#               V(s + 1bp) - V(s - 1bp)
#     Spread01 = ------------------------
#                          2
#
# The issuer curve is rebuilt after each bump so the calculation measures
# the sensitivity of the CDS to the calibrated market credit spreads.
# ============================================================================

print("\n" + LINE)
print("11. SPREAD 01 BY BUMP AND REVALUE")
print(LINE)

spread_bump = 1.0 / 10000.0

issuer_curve_up = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
    spd_bump=spread_bump,
)

issuer_curve_down = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
    spd_bump=-spread_bump,
)

value_spread_up = cds_index_contract.value(
    value_dt,
    issuer_curve_up,
    cds_recovery,
)[CLEAN]

value_spread_down = cds_index_contract.value(
    value_dt,
    issuer_curve_down,
    cds_recovery,
)[CLEAN]

spread_01_bump = (
    value_spread_up
    - value_spread_down
) / 2.0

print(
    f"{'Base Value':<35}"
    f"{clean_value:20.6f}"
)

print(
    f"{'Value +1 bp':<35}"
    f"{value_spread_up:20.6f}"
)

print(
    f"{'Value -1 bp':<35}"
    f"{value_spread_down:20.6f}"
)

print(
    f"{'Spread 01':<35}"
    f"{spread_01_bump:20.6f}"
)


# ============================================================================
# 12. CDS INDEX VALUE VERSUS RECOVERY RATE
# ============================================================================
# What this section demonstrates:
#
# Recovery affects the loss paid following default.
#
# Holding the calibrated issuer curve fixed while changing recovery isolates
# the direct valuation effect of the recovery assumption.
#
# A larger recovery rate reduces loss-given-default and therefore generally
# reduces the value of protection to a long-protection investor.
# ============================================================================

print("\n" + LINE)
print("12. CDS INDEX VALUE VERSUS RECOVERY RATE")
print(LINE)

recovery_grid = np.linspace(
    0.20,
    0.60,
    21,
)

recovery_values = []

for recovery_rate in recovery_grid:

    contract_value = cds_index_contract.value(
        value_dt,
        issuer_curve,
        recovery_rate,
    )[CLEAN]

    recovery_values.append(
        contract_value,
    )

recovery_values = np.array(
    recovery_values,
)

print(
    f"{'RECOVERY':>15}"
    f"{'CLEAN VALUE':>20}"
)

print("-" * 35)

for recovery_rate, contract_value in zip(
    recovery_grid,
    recovery_values,
):

    print(
        f"{recovery_rate:15.4f}"
        f"{contract_value:20.6f}"
    )


# ============================================================================
# 13. PLOT CDS INDEX VALUE VERSUS RECOVERY RATE
# ============================================================================

print("\n" + LINE)
print("13. CDS INDEX VALUE VERSUS RECOVERY RATE")
print(LINE)

plt.figure()

plt.plot(
    recovery_grid,
    recovery_values,
    marker="o",
)

plt.axvline(
    cds_recovery,
    linestyle="--",
    label="Base Recovery",
)

plt.axhline(
    clean_value,
    linestyle="--",
)

plt.xlabel("Recovery Rate")
plt.ylabel("Clean Value")

plt.title(
    "CDS Index Value versus Recovery Rate"
)

plt.grid(True)
plt.legend()
plt.show()


# ============================================================================
# 14. RECOVERY 01 BY BUMP AND REVALUE
# ============================================================================
# What this section demonstrates:
#
# Recovery sensitivity can also be checked directly by bumping the recovery
# assumption and revaluing the same CDS contract.
#
# Here the issuer curve is deliberately held fixed. This isolates the direct
# effect of recovery on the CDS valuation rather than recalibrating the credit
# curve after changing recovery.
# ============================================================================

print("\n" + LINE)
print("14. RECOVERY 01 BY BUMP AND REVALUE")
print(LINE)

recovery_bump = 0.01

value_recovery_up = cds_index_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery + recovery_bump,
)[CLEAN]

value_recovery_down = cds_index_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery - recovery_bump,
)[CLEAN]

recovery_01_bump = (
    value_recovery_up
    - value_recovery_down
) / 2.0

print(
    f"{'Base Recovery':<35}"
    f"{cds_recovery:20.6f}"
)

print(
    f"{'Value Recovery +1%':<35}"
    f"{value_recovery_up:20.6f}"
)

print(
    f"{'Value Recovery -1%':<35}"
    f"{value_recovery_down:20.6f}"
)

print(
    f"{'Recovery 01':<35}"
    f"{recovery_01_bump:20.6f}"
)


# ============================================================================
# 15. SUMMARY
# ============================================================================
# What this section demonstrates:
#
# Summarises the principal valuation quantities and sensitivities generated
# by the example.
# ============================================================================

print("\n" + LINE)
print("15. CDS INDEX VALUATION SUMMARY")
print(LINE)

print(
    f"{'Par Spread (bp)':<35}"
    f"{par_spread:20.6f}"
)

print(
    f"{'Contract Coupon (bp)':<35}"
    f"{index_cpn * 10000.0:20.6f}"
)

print(
    f"{'Clean Value':<35}"
    f"{clean_value:20.6f}"
)

print(
    f"{'Dirty Value':<35}"
    f"{dirty_value:20.6f}"
)

print(
    f"{'Protection Leg PV':<35}"
    f"{protection_pv:20.6f}"
)

print(
    f"{'Premium Leg PV':<35}"
    f"{premium_pv:20.6f}"
)

print(
    f"{'Spread 01':<35}"
    f"{spread_01_bump:20.6f}"
)

print(
    f"{'Recovery 01':<35}"
    f"{recovery_01_bump:20.6f}"
)

print(LINE)
