# ============================================================================
# FINANCEPY EXAMPLES - CDS Curve
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#
# This example demonstrates:
#
#   1. Construction of a CDS survival curve from market CDS spreads
#   2. Inspection of calibrated survival and default probabilities
#   3. Repricing of the CDS contracts used in the calibration
#   4. The effect of the assumed recovery rate on the calibrated curve
#
# The interest-rate curve is constructed using the common helper function
# build_ibor_curve(), keeping this example focused on credit-curve behaviour.
#
# A CDS curve combines:
#
#   - an interest-rate discount curve
#   - market CDS spreads
#   - an assumed recovery rate
#
# to infer risk-neutral survival probabilities.
#
# If Q(t) denotes survival probability to time t, then:
#
#       cumulative default probability = 1 - Q(t)
#
# Recovery affects the calibration because the protection payment following
# default depends on loss given default:
#
#       LGD = 1 - recovery rate
#
# Holding CDS spreads fixed, changing recovery therefore changes the default
# probabilities required to reproduce those spreads.
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np

from financepy.market.curves.cds_curve import CDSCurve
from financepy.products.credit.cds import CDS
from financepy.utils.calendar import (
    BusDayAdjustTypes,
    CalendarTypes,
    DateGenRuleTypes,
)
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.format_graphs import set_plot_style

from helpers import build_ibor_curve

# ============================================================================
# OUTPUT FORMAT
# ============================================================================

LINE = "=" * 100
SUBLINE = "-" * 100
set_plot_style()


# ============================================================================
# 1. CDS CURVE
# ============================================================================
#
# Construct an IBOR discount curve using the common example helper.
#
# We then create CDS contracts with maturities from one to ten years.
# Their market spreads increase from 50 bp at one year to 140 bp at ten
# years.
#
# The CDSCurve calibration determines survival probabilities consistent with
# those spreads, the discount curve and the assumed 40% recovery rate.
# ============================================================================

print("\n" + LINE)
print("1. CDS CURVE")
print(LINE)

value_dt = Date(
    20,
    12,
    2018,
)

libor_curve = build_ibor_curve(
    value_dt,
)

recovery_rate = 0.40

cds_contracts = []


print(f"{'MATURITY':>18}" f"{'SPREAD (bp)':>18}")

print(SUBLINE)


for year in range(
    1,
    11,
):

    maturity_dt = value_dt.add_months(12 * year)

    spread = 0.005 + 0.001 * (year - 1)

    cds = CDS(
        value_dt,
        maturity_dt,
        spread,
    )

    cds_contracts.append(cds)

    print(f"{str(maturity_dt):>18}" f"{spread * 10000.0:18.6f}")


issuer_curve = CDSCurve(
    value_dt,
    cds_contracts,
    libor_curve,
    recovery_rate=recovery_rate,
)


print(f"\n{'Value Date':<30}: " f"{value_dt}")

print(f"{'Recovery Rate':<30}: " f"{recovery_rate * 100.0:.2f}%")


# ============================================================================
# 2. CALIBRATED SURVIVAL PROBABILITIES
# ============================================================================
#
# CDSCurve stores the calibrated curve using:
#
#       _times      times from the valuation date
#       _qs         survival probabilities
#
# Survival probability Q(t) gives the probability, under the pricing measure,
# that the issuer survives to time t.
#
# The corresponding cumulative default probability is:
#
#       PD(t) = 1 - Q(t)
#
# Survival probabilities normally decline with horizon because there is more
# time for a default event to occur.
# ============================================================================

print("\n" + LINE)
print("2. CALIBRATED SURVIVAL PROBABILITIES")
print(LINE)

print(f"{'TIME (YEARS)':>16}" f"{'SURVIVAL PROB':>22}" f"{'DEFAULT PROB':>22}")

print(SUBLINE)


num_points = len(issuer_curve._times)


for i in range(num_points):

    time_value = issuer_curve._times[i]

    survival_probability = issuer_curve._qs[i]

    default_probability = 1.0 - survival_probability

    print(f"{time_value:16.6f}" f"{survival_probability:22.8f}" f"{default_probability:22.8f}")


# ============================================================================
# 3. REPRICE THE CALIBRATION CDS CONTRACTS
# ============================================================================
#
# Revalue the exact CDS contracts used to construct the credit curve.
#
# CDS.value() returns two values:
#
#       dirty_value, clean_value
#
# The dirty value includes accrued premium, while the clean value excludes
# accrued premium.
#
# Repricing the calibration instruments provides a useful check that the
# calibrated issuer curve is consistent with the CDS spreads used to build it.
# ============================================================================

print("\n" + LINE)
print("3. REPRICE THE CALIBRATION CDS CONTRACTS")
print(LINE)

print(f"{'CONTRACT':>10}" f"{'MATURITY':>18}" f"{'SPREAD (bp)':>18}" f"{'DIRTY VALUE':>20}" f"{'CLEAN VALUE':>20}")

print(SUBLINE)


for i, cds in enumerate(cds_contracts):

    spread = 0.005 + 0.001 * i

    dirty_value, clean_value = cds.value(
        value_dt,
        issuer_curve,
        recovery_rate,
    )

    print(
        f"{i + 1:10d}"
        f"{str(cds.maturity_dt):>18}"
        f"{spread * 10000.0:18.6f}"
        f"{dirty_value:20.8f}"
        f"{clean_value:20.8f}"
    )


# ============================================================================
# 4. PLOT THE CALIBRATED CDS CURVE
# ============================================================================
#
# Plot both survival probability and cumulative default probability.
#
# Since:
#
#       PD(t) = 1 - Q(t)
#
# the two curves contain the same information but provide complementary
# views of the calibrated credit risk.
# ============================================================================

print("\n" + LINE)
print("4. CALIBRATED CDS SURVIVAL CURVE")
print(LINE)


plot_times = np.asarray(issuer_curve._times)

plot_survival = np.asarray(issuer_curve._qs)

plot_default = 1.0 - plot_survival


plt.figure()

plt.plot(
    plot_times,
    plot_survival,
    marker="o",
    label="Survival Probability",
)

plt.plot(
    plot_times,
    plot_default,
    marker="o",
    label="Cumulative Default Probability",
)

plt.xlabel("Time (years)")

plt.ylabel("Probability")

plt.title("Calibrated CDS Survival and Default Probabilities")

plt.ylim(
    0.0,
    1.02,
)

plt.legend()
plt.grid(True)


# ============================================================================
# 5. CDS RECOVERY-RATE SENSITIVITY
# ============================================================================
#
# The second example examines how the recovery-rate assumption affects the
# calibrated survival curve.
#
# The market CDS spreads are held fixed while recovery is changed.
#
# A useful approximate relationship is:
#
#       CDS spread ~ hazard rate * (1 - recovery)
#
# Therefore:
#
#       hazard rate ~ CDS spread / (1 - recovery)
#
# For a fixed market CDS spread:
#
#   lower recovery
#       -> larger loss given default
#       -> less default intensity required
#       -> higher survival probability
#
#   higher recovery
#       -> smaller loss given default
#       -> greater default intensity required
#       -> lower survival probability
#
# FinancePy performs the full CDS calibration rather than using this simple
# approximation, but the approximation explains the direction of the effect.
# ============================================================================

print("\n" + LINE)
print("5. CDS RECOVERY-RATE SENSITIVITY")
print(LINE)


spreads = [
    0.000881720,
    0.002246440,
    0.004283100,
    0.005730380,
    0.006982450,
]

tenors = [
    "1Y",
    "3Y",
    "5Y",
    "7Y",
    "10Y",
]

effective_dt = value_dt.add_days(1)


# ============================================================================
# 5.1 CDS MARKET INPUTS
# ============================================================================

print("\nCDS market inputs")

print(f"{'TENOR':>12}" f"{'CDS SPREAD (bp)':>20}")

print(SUBLINE)


cdss = []


for tenor, spread in zip(
    tenors,
    spreads,
):

    freq_type = FrequencyTypes.MONTHLY
    accrual_dc_type = DayCountTypes.ACT_360
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.FORWARD

    cds = CDS(
        effective_dt,
        tenor,
        spread,
        notional=100,
        freq_type=freq_type,
        accrual_dc_type=accrual_dc_type,
        cal_type=cal_type,
        bd_type=bd_type,
        dg_type=dg_type,
    )

    cdss.append(cds)

    print(f"{tenor:>12}" f"{spread * 10000.0:20.6f}")


# ============================================================================
# 5.2 CALIBRATE DIFFERENT RECOVERY ASSUMPTIONS
# ============================================================================
#
# Retain every calibrated curve so that they can be compared directly.
#
# The original example successively overwrote issuer_curve. Keeping each
# curve allows us to see explicitly how recovery affects the implied
# survival probabilities.
# ============================================================================

recovery_rates = [
    0.10,
    0.575,
    0.90,
]

recovery_curves = {}


for recovery_rate in recovery_rates:

    recovery_curves[recovery_rate] = CDSCurve(
        value_dt,
        cdss,
        libor_curve,
        recovery_rate,
    )


# ============================================================================
# 6. SURVIVAL PROBABILITIES BY RECOVERY RATE
# ============================================================================

print("\n" + LINE)
print("6. SURVIVAL PROBABILITIES BY RECOVERY RATE")
print(LINE)

print(f"{'TIME':>12}" f"{'Q(t), R=10%':>20}" f"{'Q(t), R=57.5%':>20}" f"{'Q(t), R=90%':>20}")

print(SUBLINE)


reference_curve = recovery_curves[0.10]

num_points = len(reference_curve._times)


for i in range(num_points):

    time_value = reference_curve._times[i]

    q_10 = recovery_curves[0.10]._qs[i]

    q_575 = recovery_curves[0.575]._qs[i]

    q_90 = recovery_curves[0.90]._qs[i]

    print(f"{time_value:12.6f}" f"{q_10:20.8f}" f"{q_575:20.8f}" f"{q_90:20.8f}")


# ============================================================================
# 7. DEFAULT PROBABILITIES BY RECOVERY RATE
# ============================================================================
#
# Report the same calibration in terms of cumulative default probabilities:
#
#       PD(t) = 1 - Q(t)
#
# This often makes the recovery-rate effect easier to interpret.
# ============================================================================

print("\n" + LINE)
print("7. DEFAULT PROBABILITIES BY RECOVERY RATE")
print(LINE)

print(f"{'TIME':>12}" f"{'PD(t), R=10%':>20}" f"{'PD(t), R=57.5%':>20}" f"{'PD(t), R=90%':>20}")

print(SUBLINE)


for i in range(num_points):

    time_value = reference_curve._times[i]

    pd_10 = 1.0 - recovery_curves[0.10]._qs[i]

    pd_575 = 1.0 - recovery_curves[0.575]._qs[i]

    pd_90 = 1.0 - recovery_curves[0.90]._qs[i]

    print(f"{time_value:12.6f}" f"{pd_10:20.8f}" f"{pd_575:20.8f}" f"{pd_90:20.8f}")


# ============================================================================
# 8. PLOT SURVIVAL PROBABILITY VERSUS RECOVERY RATE
# ============================================================================

print("\n" + LINE)
print("8. SURVIVAL PROBABILITY VERSUS RECOVERY RATE")
print(LINE)


plt.figure()


for recovery_rate in recovery_rates:

    curve = recovery_curves[recovery_rate]

    times = np.asarray(curve._times)

    survival = np.asarray(curve._qs)

    plt.plot(
        times,
        survival,
        marker="o",
        label=(f"Recovery = " f"{recovery_rate * 100.0:.1f}%"),
    )


plt.xlabel("Time (years)")

plt.ylabel("Survival Probability")

plt.title("CDS Survival Curve versus Recovery Rate")

plt.ylim(
    0.0,
    1.02,
)

plt.legend()
plt.grid(True)


# ============================================================================
# 9. PLOT DEFAULT PROBABILITY VERSUS RECOVERY RATE
# ============================================================================

print("\n" + LINE)
print("9. DEFAULT PROBABILITY VERSUS RECOVERY RATE")
print(LINE)


plt.figure()


for recovery_rate in recovery_rates:

    curve = recovery_curves[recovery_rate]

    times = np.asarray(curve._times)

    default_probability = 1.0 - np.asarray(curve._qs)

    plt.plot(
        times,
        default_probability,
        marker="o",
        label=(f"Recovery = " f"{recovery_rate * 100.0:.1f}%"),
    )


plt.xlabel("Time (years)")

plt.ylabel("Cumulative Default Probability")

plt.title("CDS Default Probability versus Recovery Rate")

plt.ylim(
    0.0,
    1.02,
)

plt.legend()
plt.grid(True)


# ============================================================================
# 10. SUMMARY
# ============================================================================

print("\n" + LINE)
print("10. SUMMARY")
print(LINE)

print("The CDS curve converts observed market CDS spreads into " "risk-neutral survival probabilities.")

print("Survival probability generally decreases with maturity as the " "time available for default increases.")

print(
    "The CDS contracts used to calibrate the curve can be repriced "
    "against the resulting issuer curve as a calibration check."
)

print("Recovery is a calibration assumption rather than an independently " "observed quantity in this example.")

print(
    "Holding CDS spreads fixed, a higher recovery assumption implies "
    "a smaller loss given default and therefore requires greater "
    "default intensity to reproduce the same CDS spread."
)

print(
    "The recovery-rate comparison therefore produces progressively "
    "lower survival probabilities as assumed recovery increases."
)

print("\n" + LINE)
print("END OF CDS CURVE EXAMPLE")
print(LINE)


# ============================================================================
# DISPLAY ALL PLOTS
# ============================================================================

plt.show()
