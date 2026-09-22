# ============================================================================
# FINANCEPY EXAMPLES - CDS
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates:
#
#   1. Construction of an interest-rate curve
#   2. Construction of a CDS issuer curve
#   3. CDS valuation
#   4. Premium and protection legs
#   5. Clean and dirty CDS values
#   6. CDS payment schedule
#   7. Spread DV01
#   8. Interest-rate DV01
#   9. Recovery DV01
#  10. Independent bump-and-revalue verification of all three risk measures
#  11. Graphs showing the sensitivity of CDS value to:
#          - credit spreads
#          - interest rates
#          - recovery rates
#
# The important idea in the risk sections is:
#
#       Risk measure = Value(bumped market) - Value(base market)
#
# We compare this independently calculated result with the corresponding
# FinancePy risk function.
#
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_vars import CLEAN, DIRTY
from financepy.products.credit.cds import CDS
from financepy.utils.format_graphs import set_plot_style

from helpers import build_ibor_curve
from helpers import build_issuer_curve
from helpers import build_frozen_issuer_curve

# ============================================================================
# GLOBAL FORMATTING
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100
set_plot_style()


# ============================================================================
# SUPPORTING FUNCTION - BUILD INTEREST-RATE CURVE
# ============================================================================
#
# The optional bump parameter allows us to move every swap rate by the same
# amount.
#
# For example:
#
#       bump = 0.0001
#
# means that every interest rate is increased by one basis point.
#
# This is used later to independently test FinancePy's ir_dv01() function.
# ============================================================================


# ============================================================================
# SUPPORTING FUNCTION - EXTRACT CDS VALUE
# ============================================================================
#
# CDS.value() returns clean and dirty values.
#
# This helper makes the later calculations easier to read.
# ============================================================================


def dirty_value(
    cds_contract,
    value_dt,
    issuer_curve,
    recovery_rate,
):

    value = cds_contract.value(
        value_dt,
        issuer_curve,
        recovery_rate,
    )

    return value[DIRTY]


# ============================================================================
# SUPPORTING FUNCTION - SAFE PAYMENT SCHEDULE
# ============================================================================
#
# We only print payments AFTER the valuation date.
#
# This avoids asking the discount curve for discount factors at negative
# times, which would generate:
#
#       ValueError: Interpolation times must be non-negative.
#
# ============================================================================


def print_cds_payments(
    cds_contract,
    value_dt,
    issuer_curve,
):

    print(
        f"{'PAYMENT_DT':>15}" f"{'YEAR_FRAC':>14}" f"{'PAYMENT':>14}" f"{'DF':>14}" f"{'SURV_PROB':>14}" f"{'NPV':>14}"
    )

    print("-" * 85)

    num_flows = len(cds_contract.payment_dts)

    for i in range(num_flows):

        payment_dt = cds_contract.payment_dts[i]

        # Only future cash flows should be discounted from value_dt.

        if payment_dt > value_dt:

            accrual_factor = cds_contract.accrual_factors[i]
            flow = cds_contract.flows[i]

            df = issuer_curve.df(
                payment_dt,
            )

            survival_probability = issuer_curve.survival_prob(
                payment_dt,
            )

            npv = flow * df * survival_probability

            print(
                f"{str(payment_dt):>15}"
                f"{accrual_factor:14.6f}"
                f"{flow:14.6f}"
                f"{df:14.8f}"
                f"{survival_probability:14.8f}"
                f"{npv:14.6f}"
            )


# ============================================================================
# 1. MARKET SETUP
# ============================================================================

print("\n" + LINE)
print("1. MARKET SETUP")
print(LINE)

trade_dt = Date(
    15,
    8,
    2022,
)

# In this example:
#
#       value date = trade date

value_dt = trade_dt

# CDS protection normally becomes effective after the trade date.

step_in_dt = trade_dt.add_days(1)

maturity_dt = Date(
    20,
    6,
    2027,
)

cds_recovery = 0.40

notional = 1_000_000.0

cds_coupon = 0.005

long_protection = True

print(f"{'Trade Date':<40}: {trade_dt}")
print(f"{'Value Date':<40}: {value_dt}")
print(f"{'Step-In Date':<40}: {step_in_dt}")
print(f"{'Maturity Date':<40}: {maturity_dt}")

print(f"{'Notional':<40}: {notional:15,.2f}")
print(f"{'CDS Coupon':<40}: {cds_coupon * 10000.0:15.4f} bp")
print(f"{'Recovery Rate':<40}: {cds_recovery * 100.0:15.4f}%")


# ============================================================================
# 2. BUILD INTEREST-RATE CURVE
# ============================================================================

print("\n" + LINE)
print("2. BUILD INTEREST-RATE CURVE")
print(LINE)

libor_curve = build_ibor_curve(
    value_dt,
)

print("Interest-rate curve constructed.")


# ============================================================================
# 3. BUILD ISSUER CURVE
# ============================================================================

print("\n" + LINE)
print("3. BUILD ISSUER CREDIT CURVE")
print(LINE)

issuer_curve = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
)

print("Issuer CDS curve constructed.")


# ============================================================================
# 4. CREATE CDS CONTRACT
# ============================================================================

print("\n" + LINE)
print("4. CREATE CDS CONTRACT")
print(LINE)

cds_contract = CDS(
    step_in_dt,
    maturity_dt,
    cds_coupon,
    notional,
    long_protection,
)

print(f"{'Notional':<40}: {notional:15,.2f}")
print(f"{'Coupon':<40}: {cds_coupon * 10000.0:15.6f} bp")
print(f"{'Long Protection':<40}: {long_protection}")


# ============================================================================
# 5. CDS VALUATION
# ============================================================================

print("\n" + LINE)
print("5. CDS VALUATION")
print(LINE)

value = cds_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery,
)

dirty_pv = value[DIRTY]
clean_pv = value[CLEAN]

par_spread = cds_contract.par_spread(
    value_dt,
    issuer_curve,
    cds_recovery,
)

clean_price = cds_contract.clean_price(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Dirty Value':<40}: {dirty_pv:15.6f}")
print(f"{'Clean Value':<40}: {clean_pv:15.6f}")
print(f"{'Clean Price':<40}: {clean_price:15.6f}")

print(f"{'Par Spread':<40}: " f"{par_spread * 10000.0:15.6f} bp")


# ============================================================================
# 6. CDS LEGS
# ============================================================================

print("\n" + LINE)
print("6. CDS PREMIUM AND PROTECTION LEGS")
print(LINE)

protection_leg_pv = cds_contract.prot_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

premium_leg_pv = cds_contract.premium_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

dirty_rpv01, clean_rpv01 = cds_contract.rpv01(
    value_dt,
    issuer_curve,
)

print(f"{'Protection Leg PV':<40}: " f"{protection_leg_pv:15.6f}")

print(f"{'Premium Leg PV':<40}: " f"{premium_leg_pv:15.6f}")

print(f"{'Dirty RPV01':<40}: " f"{dirty_rpv01:15.8f}")

print(f"{'Clean RPV01':<40}: " f"{clean_rpv01:15.8f}")


# ============================================================================
# 7. ACCRUED PREMIUM
# ============================================================================

print("\n" + LINE)
print("7. CDS ACCRUED PREMIUM")
print(LINE)

accrued_days = cds_contract.accrued_days(
    value_dt,
)

accrued_interest = cds_contract.accrued_interest(
    value_dt,
)

print(f"{'Accrued Days':<40}: " f"{accrued_days}")

print(f"{'Accrued Premium':<40}: " f"{accrued_interest:15.6f}")


# ============================================================================
# 8. PAYMENT SCHEDULE
# ============================================================================
#
# Only future payments are shown.
#
# For each payment we display:
#
#       accrual fraction
#       contractual premium payment
#       discount factor
#       survival probability
#       discounted survival-weighted payment
#
# ============================================================================

print("\n" + LINE)
print("8. CDS PAYMENT SCHEDULE")
print(LINE)

print_cds_payments(
    cds_contract,
    value_dt,
    issuer_curve,
)


# ============================================================================
# 9. FINANCEPY RISK MEASURES
# ============================================================================

print("\n" + LINE)
print("9. FINANCEPY CDS RISK MEASURES")
print(LINE)

spread_dv01_function = cds_contract.spread_dv01(
    value_dt,
    issuer_curve,
    cds_recovery,
)

ir_dv01_function = cds_contract.ir_dv01(
    value_dt,
    issuer_curve,
    cds_recovery,
)

recovery_dv01_function = cds_contract.recovery_dv01(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Spread DV01':<40}: " f"{spread_dv01_function:15.8f}")

print(f"{'IR DV01':<40}: " f"{ir_dv01_function:15.8f}")

print(f"{'Recovery DV01':<40}: " f"{recovery_dv01_function:15.8f}")


# ============================================================================
# 10. TEST SPREAD DV01 BY BUMP AND REVALUE
# ============================================================================
#
# Spread DV01 measures the change in CDS value caused by a one basis point
# increase in the market CDS spreads used to calibrate the issuer curve.
#
# We independently reproduce it:
#
#       1. Start with the base CDS value.
#       2. Increase every calibration CDS spread by 1 bp.
#       3. Rebuild the issuer curve.
#       4. Revalue the CDS.
#       5. Calculate:
#
#              bumped value - base value
#
# ============================================================================

print("\n" + LINE)
print("10. TEST SPREAD DV01 BY BUMP AND REVALUE")
print(LINE)

spread_bump = 1.0 / 10000.0

base_value = dirty_value(
    cds_contract,
    value_dt,
    issuer_curve,
    cds_recovery,
)

issuer_curve_spread_up = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
    spd_bump=spread_bump,
)

spread_bumped_value = dirty_value(
    cds_contract,
    value_dt,
    issuer_curve_spread_up,
    cds_recovery,
)

spread_dv01_bump = spread_bumped_value - base_value

spread_error = spread_dv01_function - spread_dv01_bump

print(f"{'Base Value':<45}: {base_value:15.8f}")

print(f"{'Value after +1 bp CDS spread bump':<45}: " f"{spread_bumped_value:15.8f}")

print(f"{'Bump-and-Revalue Spread DV01':<45}: " f"{spread_dv01_bump:15.8f}")

print(f"{'FinancePy spread_dv01()':<45}: " f"{spread_dv01_function:15.8f}")

print(f"{'Difference':<45}: " f"{spread_error:15.10f}")


# ============================================================================
# 11. TEST IR DV01 BY BUMP AND REVALUE
# ============================================================================
#
# Interest-rate DV01 measures the change in CDS value following a parallel
# one basis point increase in the interest-rate curve.
#
# There are actually TWO mechanisms through which the rate bump can affect
# the CDS value.
#
#
# EFFECT 1 - DIRECT DISCOUNTING
# -----------------------------
#
# Changing interest rates changes the discount factors:
#
#                       DF(t)
#
# used to present-value both:
#
#       - the premium leg
#       - the protection leg
#
# To isolate this effect we:
#
#       1. bump the interest-rate curve
#       2. keep the survival probabilities unchanged
#       3. revalue the CDS
#
#
# EFFECT 2 - ISSUER CURVE RECALIBRATION
# -------------------------------------
#
# Market CDS spreads are calibration instruments.
#
# Once the discount curve changes, the original survival probabilities will
# generally no longer reproduce exactly the same market CDS spreads.
#
# Therefore the issuer curve is recalibrated using the bumped discount curve.
#
#
# We therefore calculate three values:
#
#
#       V0 = PV(base rates, base survival curve)
#
#       V1 = PV(bumped rates, frozen survival curve)
#
#       V2 = PV(bumped rates, recalibrated survival curve)
#
#
# This gives:
#
#
#       Direct Discounting Effect
#
#           = V1 - V0
#
#
#       Credit Recalibration Effect
#
#           = V2 - V1
#
#
#       Total IR DV01
#
#           = V2 - V0
#
#
# and therefore:
#
#
#       Total IR DV01
#
#           = Direct Discounting Effect
#             + Credit Recalibration Effect
#
#
# The final total is compared with FinancePy's ir_dv01().
# ============================================================================

print("\n" + LINE)
print("11. TEST IR DV01 BY BUMP AND REVALUE")
print(LINE)


# ============================================================================
# 11.1 BASE VALUE
# ============================================================================

print("\n" + SUBLINE)
print("11.1 BASE CDS VALUE")
print(SUBLINE)

ir_bump = 1.0 / 10000.0

v0 = dirty_value(
    cds_contract,
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(
    f"{'V0 - Base CDS Value':<50}: "
    f"{v0:15.8f}"
)


# ============================================================================
# 11.2 BUMP THE INTEREST-RATE CURVE
# ============================================================================

print("\n" + SUBLINE)
print("11.2 BUMP INTEREST-RATE CURVE BY +1 BP")
print(SUBLINE)

libor_curve_ir_up = build_ibor_curve(
    value_dt,
    bump=ir_bump,
)

print(
    f"{'Parallel Interest-Rate Bump':<50}: "
    f"{ir_bump * 10000.0:15.6f} bp"
)


# ============================================================================
# 11.3 DIRECT DISCOUNTING EFFECT
# ============================================================================
#
# Replace the discount curve but DO NOT recalibrate survival probabilities.
#
# This isolates the effect of changing discount factors.
# ============================================================================

print("\n" + SUBLINE)
print("11.3 DIRECT DISCOUNTING EFFECT")
print(SUBLINE)

frozen_issuer_curve = build_frozen_issuer_curve(
    issuer_curve,
    libor_curve_ir_up,
    value_dt,
    cds_recovery,
)

v1 = dirty_value(
    cds_contract,
    value_dt,
    frozen_issuer_curve,
    cds_recovery,
)

direct_discounting_effect = (
    v1
    - v0
)

print(
    f"{'V0 - Base Value':<50}: "
    f"{v0:15.8f}"
)

print(
    f"{'V1 - Bumped Rates / Frozen Credit':<50}: "
    f"{v1:15.8f}"
)

print(
    f"{'Direct Discounting Effect (V1 - V0)':<50}: "
    f"{direct_discounting_effect:15.8f}"
)


# ============================================================================
# 11.4 RECALIBRATE THE ISSUER CURVE
# ============================================================================
#
# We now rebuild the issuer curve using:
#
#       - the same market CDS spreads
#       - the bumped interest-rate curve
#
# FinancePy must adjust the survival probabilities so that the calibration
# CDS instruments once again reproduce their market spreads.
# ============================================================================

print("\n" + SUBLINE)
print("11.4 ISSUER CURVE RECALIBRATION EFFECT")
print(SUBLINE)

issuer_curve_ir_up = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve_ir_up,
    cds_recovery,
)

v2 = dirty_value(
    cds_contract,
    value_dt,
    issuer_curve_ir_up,
    cds_recovery,
)

credit_recalibration_effect = (
    v2
    - v1
)

print(
    f"{'V1 - Bumped Rates / Frozen Credit':<50}: "
    f"{v1:15.8f}"
)

print(
    f"{'V2 - Bumped Rates / Recalibrated Credit':<50}: "
    f"{v2:15.8f}"
)

print(
    f"{'Credit Recalibration Effect (V2 - V1)':<50}: "
    f"{credit_recalibration_effect:15.8f}"
)


# ============================================================================
# 11.5 TOTAL BUMP-AND-REVALUE IR DV01
# ============================================================================

print("\n" + SUBLINE)
print("11.5 TOTAL IR DV01")
print(SUBLINE)

ir_dv01_bump = (
    v2
    - v0
)

ir_effect_sum = (
    direct_discounting_effect
    + credit_recalibration_effect
)

print(
    f"{'Direct Discounting Effect':<50}: "
    f"{direct_discounting_effect:15.8f}"
)

print(
    f"{'Credit Recalibration Effect':<50}: "
    f"{credit_recalibration_effect:15.8f}"
)

print(SUBLINE)

print(
    f"{'Sum of Two Effects':<50}: "
    f"{ir_effect_sum:15.8f}"
)

print(
    f"{'Total Bump-and-Revalue IR DV01 (V2 - V0)':<50}: "
    f"{ir_dv01_bump:15.8f}"
)


# ============================================================================
# 11.6 COMPARE WITH FINANCEPY IR DV01
# ============================================================================

print("\n" + SUBLINE)
print("11.6 COMPARE WITH FINANCEPY ir_dv01()")
print(SUBLINE)

ir_dv01_function = cds_contract.ir_dv01(
    value_dt,
    issuer_curve,
    cds_recovery,
)

ir_error = (
    ir_dv01_function
    - ir_dv01_bump
)

print(
    f"{'FinancePy ir_dv01()':<50}: "
    f"{ir_dv01_function:15.8f}"
)

print(
    f"{'Manual Bump-and-Revalue IR DV01':<50}: "
    f"{ir_dv01_bump:15.8f}"
)

print(
    f"{'Difference':<50}: "
    f"{ir_error:15.10f}"
)

ir_test = np.isclose(
    ir_dv01_function,
    ir_dv01_bump,
    rtol=1.0e-5,
    atol=1.0e-6,
)

print(
    f"{'IR DV01 TEST':<50}: "
    f"{'PASS' if ir_test else 'FAIL'}"
)


# ============================================================================
# 11.7 CHECK THE DECOMPOSITION
# ============================================================================

print("\n" + SUBLINE)
print("11.7 CHECK IR DV01 DECOMPOSITION")
print(SUBLINE)

decomposition_error = (
    ir_dv01_bump
    - ir_effect_sum
)

print(
    f"{'Total IR DV01':<50}: "
    f"{ir_dv01_bump:15.8f}"
)

print(
    f"{'Direct + Recalibration':<50}: "
    f"{ir_effect_sum:15.8f}"
)

print(
    f"{'Difference':<50}: "
    f"{decomposition_error:15.10f}"
)

decomposition_test = np.isclose(
    ir_dv01_bump,
    ir_effect_sum,
    rtol=1.0e-12,
    atol=1.0e-12,
)

print(
    f"{'DECOMPOSITION TEST':<50}: "
    f"{'PASS' if decomposition_test else 'FAIL'}"
)


# ============================================================================
# 11.8 SHOW HOW THE SURVIVAL CURVE CHANGES
# ============================================================================
#
# Compare the original calibrated survival probabilities with those obtained
# after the +1 bp interest-rate bump.
#
# The market CDS spreads themselves have NOT changed.
#
# Any change in Q(t) therefore comes from recalibration to the new discount
# curve.
# ============================================================================

print("\n" + SUBLINE)
print("11.8 SURVIVAL PROBABILITIES BEFORE AND AFTER RATE BUMP")
print(SUBLINE)

base_times = np.asarray(
    issuer_curve._times,
)

base_qs = np.asarray(
    issuer_curve._qs,
)

bumped_times = np.asarray(
    issuer_curve_ir_up._times,
)

bumped_qs = np.asarray(
    issuer_curve_ir_up._qs,
)

print(
    f"{'TIME':>12}"
    f"{'BASE Q(t)':>18}"
    f"{'BUMPED Q(t)':>18}"
    f"{'CHANGE':>18}"
)

print("-" * 66)

for (
    time,
    base_q,
    bumped_q,
) in zip(
    base_times,
    base_qs,
    bumped_qs,
):

    q_change = (
        bumped_q
        - base_q
    )

    print(
        f"{time:12.6f}"
        f"{base_q:18.10f}"
        f"{bumped_q:18.10f}"
        f"{q_change:18.10f}"
    )


# ============================================================================
# 11.9 GRAPH - IR DV01 DECOMPOSITION
# ============================================================================

print("\n" + SUBLINE)
print("11.9 PLOT IR DV01 DECOMPOSITION")
print(SUBLINE)

effect_names = [
    "Direct\nDiscounting",
    "Credit Curve\nRecalibration",
    "Total\nIR DV01",
]

effect_values = [
    direct_discounting_effect,
    credit_recalibration_effect,
    ir_dv01_bump,
]

plt.figure(
    figsize=(9, 6),
)

plt.bar(
    effect_names,
    effect_values,
)

plt.axhline(
    0.0,
    linestyle="--",
)

plt.ylabel(
    "Change in CDS Value"
)

plt.title(
    "Decomposition of CDS Interest-Rate DV01"
)

plt.grid(
    True,
    axis="y",
)

plt.tight_layout()
plt.show()


# ============================================================================
# 11.10 GRAPH - SURVIVAL CURVE BEFORE AND AFTER RATE BUMP
# ============================================================================

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    base_times,
    base_qs,
    marker="o",
    label="Base Interest Rates",
)

plt.plot(
    bumped_times,
    bumped_qs,
    marker="o",
    label="Interest Rates +1 bp",
)

plt.xlabel(
    "Time (years)"
)

plt.ylabel(
    "Survival Probability"
)

plt.title(
    "Effect of Interest-Rate Bump on Calibrated CDS Survival Curve"
)

plt.ylim(
    0.0,
    1.02,
)

plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ============================================================================
# 11.11 INTERPRETATION
# ============================================================================

print("\n" + SUBLINE)
print("11.11 INTERPRETATION")
print(SUBLINE)

print(
    """
The interest-rate sensitivity can now be understood as two separate effects.

V0
--
The original CDS value using:

    base interest-rate curve
    base calibrated survival curve


V1
--
The CDS value after increasing interest rates by one basis point while
holding the survival probabilities fixed.

Therefore:

    V1 - V0

is the DIRECT DISCOUNTING EFFECT.


V2
--
The CDS value after increasing interest rates by one basis point and then
recalibrating the issuer survival curve to the original market CDS spreads.

Therefore:

    V2 - V1

is the ISSUER CURVE RECALIBRATION EFFECT.


Finally:

    V2 - V0

is the complete bump-and-revalue IR DV01.

By construction:

    V2 - V0

        = (V1 - V0) + (V2 - V1)

so:

    TOTAL IR DV01

        = DIRECT DISCOUNTING EFFECT

        + ISSUER RECALIBRATION EFFECT


The FinancePy ir_dv01() result should agree with V2 - V0 if the manual bump
uses the same bump size and recalibration convention as the library function.
"""
)

# ============================================================================
# 12. TEST RECOVERY DV01 BY BUMP AND REVALUE
# ============================================================================
#
# Recovery DV01 measures the effect of changing the recovery assumption.
#
# Changing recovery changes:
#
#           Loss Given Default = 1 - Recovery
#
# and therefore changes the hazard rates implied by the CDS market spreads.
#
# Consequently the issuer curve MUST be rebuilt after changing recovery.
#
# ============================================================================

print("\n" + LINE)
print("12. TEST RECOVERY DV01 BY BUMP AND REVALUE")
print(LINE)

recovery_bump = 0.01

recovery_rate_up = cds_recovery + recovery_bump

issuer_curve_recovery_up = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    recovery_rate_up,
)

recovery_bumped_value = dirty_value(
    cds_contract,
    value_dt,
    issuer_curve_recovery_up,
    recovery_rate_up,
)

recovery_dv01_bump = recovery_bumped_value - base_value

recovery_error = recovery_dv01_function - recovery_dv01_bump

print(f"{'Base Recovery Rate':<45}: " f"{cds_recovery * 100.0:14.6f}%")

print(f"{'Bumped Recovery Rate':<45}: " f"{recovery_rate_up * 100.0:14.6f}%")

print(f"{'Base Value':<45}: " f"{base_value:15.8f}")

print(f"{'Value after recovery bump':<45}: " f"{recovery_bumped_value:15.8f}")

print(f"{'Bump-and-Revalue Recovery DV01':<45}: " f"{recovery_dv01_bump:15.8f}")

print(f"{'FinancePy recovery_dv01()':<45}: " f"{recovery_dv01_function:15.8f}")

print(f"{'Difference':<45}: " f"{recovery_error:15.10f}")


# ============================================================================
# 13. RISK TEST SUMMARY
# ============================================================================

print("\n" + LINE)
print("13. RISK MEASURE TEST SUMMARY")
print(LINE)

print(f"{'RISK':<20}" f"{'FINANCEPY':>20}" f"{'BUMP/REVALUE':>20}" f"{'DIFFERENCE':>20}")

print(SUBLINE)

risk_tests = [
    (
        "Spread DV01",
        spread_dv01_function,
        spread_dv01_bump,
    ),
    (
        "IR DV01",
        ir_dv01_function,
        ir_dv01_bump,
    ),
    (
        "Recovery DV01",
        recovery_dv01_function,
        recovery_dv01_bump,
    ),
]

for (
    risk_name,
    function_value,
    bump_value,
) in risk_tests:

    difference = function_value - bump_value

    print(f"{risk_name:<20}" f"{function_value:20.8f}" f"{bump_value:20.8f}" f"{difference:20.10f}")


# ============================================================================
# 14. ASSERTION TESTS
# ============================================================================
#
# The FinancePy functions and our manual bump calculations should agree.
#
# We use np.isclose rather than exact equality because these calculations
# involve numerical curve calibration and interpolation.
# ============================================================================

print("\n" + LINE)
print("14. AUTOMATED TESTS")
print(LINE)

spread_test = np.isclose(
    spread_dv01_function,
    spread_dv01_bump,
    rtol=1.0e-5,
    atol=1.0e-6,
)

ir_test = np.isclose(
    ir_dv01_function,
    ir_dv01_bump,
    rtol=1.0e-5,
    atol=1.0e-6,
)

recovery_test = np.isclose(
    recovery_dv01_function,
    recovery_dv01_bump,
    rtol=1.0e-5,
    atol=1.0e-6,
)

print(f"{'Spread DV01 test':<40}: " f"{'PASS' if spread_test else 'FAIL'}")

print(f"{'IR DV01 test':<40}: " f"{'PASS' if ir_test else 'FAIL'}")

print(f"{'Recovery DV01 test':<40}: " f"{'PASS' if recovery_test else 'FAIL'}")


# ============================================================================
# 15. CDS VALUE VERSUS ISSUER SPREAD
# ============================================================================
#
# We now go beyond a single 1 bp bump and examine a range of parallel spread
# movements.
#
# For a buyer of protection, increasing ISSUER spreads should generally make
# an existing fixed-coupon protection position more valuable.
#
# ============================================================================

print("\n" + LINE)
print("15. CDS VALUE VERSUS ISSUER SPREAD")
print(LINE)

spread_bumps_bp = np.arange(
    -5,
    51,
    5,
)

spread_values = []

for bump_bp in spread_bumps_bp:

    bump = bump_bp / 10000.0

    bumped_curve = build_issuer_curve(
        value_dt,
        step_in_dt,
        libor_curve,
        cds_recovery,
        spd_bump=bump,
    )

    bumped_value = dirty_value(
        cds_contract,
        value_dt,
        bumped_curve,
        cds_recovery,
    )

    spread_values.append(
        bumped_value,
    )

spread_values = np.asarray(
    spread_values,
)

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    spread_bumps_bp,
    spread_values,
    marker="o",
)

plt.axvline(
    0.0,
    linestyle="--",
)

plt.axhline(
    base_value,
    linestyle="--",
)

plt.xlabel(
    "Parallel CDS Spread Bump (bp)",
)

plt.ylabel(
    "CDS Dirty Value",
)

plt.title("CDS Value versus ISSUER Spread")

plt.grid(True)
plt.show()


# ============================================================================
# 16. CDS VALUE VERSUS INTEREST RATES
# ============================================================================
#
# Here every swap rate used to construct the discount curve is shifted by the
# same amount.
#
# For every bumped interest-rate curve we recalibrate the issuer CDS curve.
#
# ============================================================================

print("\n" + LINE)
print("16. CDS VALUE VERSUS INTEREST RATES")
print(LINE)

rate_bumps_bp = np.arange(
    -5,
    51,
    5,
)

rate_values = []

for bump_bp in rate_bumps_bp:

    bump = bump_bp / 10000.0

    bumped_libor_curve = build_ibor_curve(
        value_dt,
        bump=bump,
    )

    bumped_issuer_curve = build_issuer_curve(
        value_dt,
        step_in_dt,
        bumped_libor_curve,
        cds_recovery
    )

    bumped_value = dirty_value(
        cds_contract,
        value_dt,
        bumped_issuer_curve,
        cds_recovery
    )

    rate_values.append(
        bumped_value,
    )

rate_values = np.asarray(
    rate_values,
)

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    rate_bumps_bp,
    rate_values,
    marker="o",
)

plt.axvline(
    0.0,
    linestyle="--",
)

plt.axhline(
    base_value,
    linestyle="--",
)

plt.xlabel(
    "Parallel Interest-Rate Bump (bp)",
)

plt.ylabel(
    "CDS Dirty Value",
)

plt.title("CDS Value versus Interest Rates")

plt.grid(True)
plt.show()


# ============================================================================
# 17. CDS VALUE VERSUS RECOVERY RATE
# ============================================================================
#
# Finally we vary the recovery assumption.
#
# At every recovery rate we rebuild the CDS issuer curve so that the original
# market CDS spreads remain the calibration instruments.
#
# This is important because changing recovery without recalibrating the
# hazard curve would represent a different risk experiment.
#
# ============================================================================

print("\n" + LINE)
print("17. CDS VALUE VERSUS RECOVERY RATE")
print(LINE)

recovery_rates = np.linspace(
    0.10,
    0.70,
    25,
)

recovery_values = []

for recovery_rate in recovery_rates:

    bumped_issuer_curve = build_issuer_curve(
        value_dt,
        step_in_dt,
        libor_curve,
        recovery_rate,
    )

    bumped_value = dirty_value(
        cds_contract,
        value_dt,
        bumped_issuer_curve,
        recovery_rate,
    )

    recovery_values.append(
        bumped_value,
    )

recovery_values = np.asarray(
    recovery_values,
)

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    recovery_rates * 100.0,
    recovery_values,
    marker="o",
)

plt.axvline(
    cds_recovery * 100.0,
    linestyle="--",
)

plt.axhline(
    base_value,
    linestyle="--",
)

plt.xlabel(
    "Recovery Rate (%)",
)

plt.ylabel(
    "CDS Dirty Value",
)

plt.title("CDS Value versus Recovery Rate")

plt.grid(True)
plt.show()


# ============================================================================
# 18. SURVIVAL PROBABILITY CURVE
# ============================================================================
#
# The CDS calibration converts market CDS spreads into a term structure of
# survival probabilities.
#
# Q(t) is the risk-neutral probability of surviving from today until time t.
#
# Therefore:
#
#       Q(0) = 1
#
# and Q(t) normally decreases as the horizon increases.
#
# ============================================================================

print("\n" + LINE)
print("18. CALIBRATED SURVIVAL PROBABILITY CURVE")
print(LINE)

survival_times = np.asarray(
    issuer_curve._times,
)

survival_probabilities = np.asarray(
    issuer_curve._qs,
)

print(f"{'TIME':>15}" f"{'SURVIVAL PROBABILITY':>25}")

print("-" * 40)

for (
    time,
    survival_probability,
) in zip(
    survival_times,
    survival_probabilities,
):

    print(f"{time:15.6f}" f"{survival_probability:25.8f}")

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    survival_times,
    survival_probabilities,
    marker="o",
)

plt.xlabel(
    "Time (years)",
)

plt.ylabel(
    "Survival Probability",
)

plt.title("Calibrated CDS Survival Curve")

plt.ylim(
    0.0,
    1.02,
)

plt.grid(True)
plt.show()


# ============================================================================
# 19. INTERPRETATION
# ============================================================================

print("\n" + LINE)
print("19. INTERPRETATION")
print(LINE)

print("""
SPREAD DV01
-----------
Spread DV01 measures the change in CDS value caused by a one basis point
parallel increase in the CDS spreads used to calibrate the issuer curve.

The independent test is:

    spread DV01
        = PV(spreads + 1 bp)
        - PV(base spreads)


INTEREST-RATE DV01
------------------
IR DV01 measures the change in CDS value caused by a one basis point
parallel increase in the interest-rate curve.

Because the issuer hazard curve is calibrated using discount factors, the
issuer curve must be rebuilt after the interest-rate curve is bumped.

The independent test is:

    IR DV01
        = PV(rates + 1 bp, recalibrated issuer curve)
        - PV(base rates)


RECOVERY DV01
-------------
Recovery affects loss given default:

    LGD = 1 - Recovery

However, market CDS spreads are held fixed in the recovery sensitivity
calculation. Consequently the issuer hazard curve must be recalibrated after
the recovery assumption changes.

The independent test is:

    Recovery DV01
        = PV(bumped recovery, recalibrated issuer curve)
        - PV(base recovery)


WHY THESE TESTS ARE USEFUL
--------------------------
The built-in FinancePy risk functions and the manual bump-and-revalue
calculations are obtained through separate calculation paths.

Agreement therefore provides a useful numerical check that:

    1. the market input being bumped is the intended risk factor,
    2. the appropriate curves are being recalibrated,
    3. the CDS is being revalued consistently,
    4. the sign and units of the reported sensitivity are understood.
""")
