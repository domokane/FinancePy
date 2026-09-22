# ============================================================================
# FINANCEPY EXAMPLES - CDSOption
# ============================================================================
#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import numpy as np
import matplotlib.pyplot as plt

from financepy.utils.global_vars import CLEAN, DIRTY
from financepy.utils.date import Date
from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_option import CDSOption
from financepy.utils.format_graphs import set_plot_style

from helpers import build_ibor_curve
from helpers import build_issuer_curve

# ============================================================================
# 1. CDS OPTION
# ============================================================================
#
# What this section demonstrates:
#
# - Builds the interest-rate discount curve.
# - Builds the issuer CDS survival curve.
# - Values a spot CDS.
# - Reports clean and dirty CDS values.
# - Reports protection and premium leg PVs.
# - Calculates RPV01.
# - Constructs a forward-starting CDS.
# - Prices CDS options across a range of strikes.
# - Recovers implied volatility from each calculated option value.
#
# Date convention:
#
#     trade_dt   = valuation date
#     value_dt   = trade date
#     step_in_dt = trade date + 1 day
#
# ============================================================================

print("\n" + "=" * 78)
print("1. CDS OPTION")
print("=" * 78)
set_plot_style()


# ============================================================================
# 1.1 DATES
# ============================================================================

trade_dt = Date(
    5,
    2,
    2014,
)

value_dt = trade_dt

step_in_dt = trade_dt.add_days(
    1,
)

expiry_dt = Date(
    20,
    3,
    2014,
)

maturity_dt = Date(
    20,
    6,
    2019,
)


# ============================================================================
# 1.2 BUILD INTEREST-RATE CURVE
# ============================================================================

libor_curve = build_ibor_curve(
    value_dt,
)


# ============================================================================
# 1.3 BUILD ISSUER CREDIT CURVE
# ============================================================================

cds_recovery = 0.40

issuer_curve = build_issuer_curve(
    value_dt,
    step_in_dt,
    libor_curve,
    cds_recovery,
)


# ============================================================================
# 2. SPOT CDS
# ============================================================================
#
# First examine the underlying CDS used by the option framework.
#
# The CDS coupon is set to zero here so that the protection-leg economics can
# be inspected directly.
# ============================================================================

print("\n" + "=" * 78)
print("2. SPOT CDS")
print("=" * 78)

notional = 100.0

long_protection = False

cds_cpn = 0.0

cds_contract = CDS(
    step_in_dt,
    maturity_dt,
    cds_cpn,
    notional,
    long_protection,
)


# ============================================================================
# 2.1 PAR SPREAD
# ============================================================================

spd = cds_contract.par_spread(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Par Spread':<35}: " f"{spd * 10000.0:14.6f} bp")


# ============================================================================
# 2.2 CLEAN AND DIRTY VALUE
# ============================================================================

value = cds_contract.value(
    value_dt,
    issuer_curve,
    cds_recovery,
)

dirty_value = value[DIRTY]
clean_value = value[CLEAN]

print(f"{'Dirty Value':<35}: " f"{dirty_value:14.8f}")

print(f"{'Clean Value':<35}: " f"{clean_value:14.8f}")


# ============================================================================
# 2.3 CLEAN PRICE
# ============================================================================

clean_price = cds_contract.clean_price(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Clean Price':<35}: " f"{clean_price:14.8f}")


# ============================================================================
# 2.4 ACCRUED COUPON
# ============================================================================

accrued_days = cds_contract.accrued_days(
    value_dt,
)

accrued_interest = cds_contract.accrued_interest(
    value_dt,
)

print(f"{'Accrued Days':<35}: " f"{accrued_days}")

print(f"{'Accrued Coupon':<35}: " f"{accrued_interest:14.8f}")


# ============================================================================
# 2.5 PROTECTION LEG
# ============================================================================

prot_pv = cds_contract.prot_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Protection Leg PV':<35}: " f"{prot_pv:14.8f}")


# ============================================================================
# 2.6 PREMIUM LEG
# ============================================================================

prem_pv = cds_contract.premium_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Premium Leg PV':<35}: " f"{prem_pv:14.8f}")


# ============================================================================
# 2.7 RPV01
# ============================================================================

dirty_rpv01, clean_rpv01 = cds_contract.rpv01(
    value_dt,
    issuer_curve,
)

print(f"{'Dirty RPV01':<35}: " f"{dirty_rpv01:14.8f}")

print(f"{'Clean RPV01':<35}: " f"{clean_rpv01:14.8f}")


# ============================================================================
# 3. FORWARD CDS
# ============================================================================
#
# The CDS underlying the option begins at the option expiry date.
#
# This provides the forward CDS spread and valuation quantities relevant to
# the CDS option.
# ============================================================================

print("\n" + "=" * 78)
print("3. FORWARD CDS")
print("=" * 78)

forward_cds = CDS(
    expiry_dt,
    maturity_dt,
    cds_cpn,
    notional,
    long_protection,
)


# ============================================================================
# 3.1 FORWARD PAR SPREAD
# ============================================================================

forward_spd = forward_cds.par_spread(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Forward Par Spread':<35}: " f"{forward_spd * 10000.0:14.6f} bp")


# ============================================================================
# 3.2 FORWARD CLEAN AND DIRTY VALUE
# ============================================================================

forward_value = forward_cds.value(
    value_dt,
    issuer_curve,
    cds_recovery,
)

forward_dirty_value = forward_value[DIRTY]
forward_clean_value = forward_value[CLEAN]

print(f"{'Forward Dirty Value':<35}: " f"{forward_dirty_value:14.8f}")

print(f"{'Forward Clean Value':<35}: " f"{forward_clean_value:14.8f}")


# ============================================================================
# 3.3 FORWARD PROTECTION LEG
# ============================================================================

forward_prot_pv = forward_cds.prot_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Forward Protection Leg PV':<35}: " f"{forward_prot_pv:14.8f}")


# ============================================================================
# 3.4 FORWARD PREMIUM LEG
# ============================================================================

forward_prem_pv = forward_cds.premium_leg_pv(
    value_dt,
    issuer_curve,
    cds_recovery,
)

print(f"{'Forward Premium Leg PV':<35}: " f"{forward_prem_pv:14.8f}")


# ============================================================================
# 3.5 FORWARD RPV01
# ============================================================================

forward_dirty_rpv01, forward_clean_rpv01 = forward_cds.rpv01(
    value_dt,
    issuer_curve,
)

print(f"{'Forward Dirty RPV01':<35}: " f"{forward_dirty_rpv01:14.8f}")

print(f"{'Forward Clean RPV01':<35}: " f"{forward_clean_rpv01:14.8f}")


# ============================================================================
# 4. CDS OPTIONS
# ============================================================================
#
# Price CDS options over a range of strikes.
#
# Long protection:
#
#     benefits when credit spreads widen.
#
# Short protection:
#
#     benefits when credit spreads tighten.
#
# The same volatility is used to generate every option price. The implied
# volatility calculation then solves backwards for the volatility that
# reproduces each calculated price.
# ============================================================================

print("\n" + "=" * 78)
print("4. CDS OPTIONS")
print("=" * 78)

cds_cpn = 0.01

volatility = 0.30

print(f"{'Value Date':<35}: " f"{str(value_dt)}")

print(f"{'Step-In Date':<35}: " f"{str(step_in_dt)}")

print(f"{'Expiry Date':<35}: " f"{str(expiry_dt)}")

print(f"{'Maturity Date':<35}: " f"{str(maturity_dt)}")

print(f"{'CDS Coupon':<35}: " f"{cds_cpn * 100.0:13.6f}%")

print(f"{'Input Volatility':<35}: " f"{volatility * 100.0:13.6f}%")


# ============================================================================
# 4.1 STRIKE GRID
# ============================================================================

strikes = np.linspace(
    100.0,
    300.0,
    41,
)

long_values = []
long_implied_vols = []

short_values = []
short_implied_vols = []


# ============================================================================
# 4.2 LONG-PROTECTION OPTIONS
# ============================================================================

print("\n" + "-" * 78)
print("LONG-PROTECTION CDS OPTIONS")
print("-" * 78)

print(f"{'STRIKE (bp)':>14}" f"{'OPTION VALUE':>20}" f"{'IMPLIED VOL (%)':>22}")

print("-" * 56)

for strike in strikes:

    long_protection = True

    cds_option = CDSOption(
        expiry_dt,
        maturity_dt,
        strike / 10000.0,
        notional,
        long_protection,
    )

    option_value = cds_option.value(
        value_dt,
        issuer_curve,
        volatility,
    )

    implied_vol = cds_option.implied_volatility(
        value_dt,
        issuer_curve,
        option_value,
    )

    long_values.append(
        option_value,
    )

    long_implied_vols.append(
        implied_vol,
    )

    print(f"{strike:14.4f}" f"{option_value:20.8f}" f"{implied_vol * 100.0:22.8f}")

print("-" * 56)


# ============================================================================
# 4.3 SHORT-PROTECTION OPTIONS
# ============================================================================

print("\n" + "-" * 78)
print("SHORT-PROTECTION CDS OPTIONS")
print("-" * 78)

print(f"{'STRIKE (bp)':>14}" f"{'OPTION VALUE':>20}" f"{'IMPLIED VOL (%)':>22}")

print("-" * 56)

for strike in strikes:

    long_protection = False

    cds_option = CDSOption(
        expiry_dt,
        maturity_dt,
        strike / 10000.0,
        notional,
        long_protection,
    )

    option_value = cds_option.value(
        value_dt,
        issuer_curve,
        volatility,
    )

    implied_vol = cds_option.implied_volatility(
        value_dt,
        issuer_curve,
        option_value,
    )

    short_values.append(
        option_value,
    )

    short_implied_vols.append(
        implied_vol,
    )

    print(f"{strike:14.4f}" f"{option_value:20.8f}" f"{implied_vol * 100.0:22.8f}")

print("-" * 56)


# ============================================================================
# 5. CDS OPTION VALUE VERSUS STRIKE
# ============================================================================
#
# The strike is the spread at which the option holder can enter the underlying
# CDS.
#
# For a long-protection option, increasing the strike makes the right to buy
# protection less attractive and therefore reduces option value.
#
# For a short-protection option, the opposite strike dependence is expected.
#
# Plotting both curves shows the asymmetric exposure of protection buyers and
# protection sellers to the CDS spread.
# ============================================================================

print("\n" + "=" * 78)
print("5. CDS OPTION VALUE VERSUS STRIKE")
print("=" * 78)

long_values = np.asarray(
    long_values,
)

short_values = np.asarray(
    short_values,
)

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    strikes,
    long_values,
    marker="o",
    label="Long Protection",
)

plt.plot(
    strikes,
    short_values,
    marker="o",
    label="Short Protection",
)

plt.axvline(
    forward_spd * 10000.0,
    linestyle="--",
    label="Forward CDS Spread",
)

plt.xlabel(
    "Strike Spread (bp)",
)

plt.ylabel(
    "Option Value",
)

plt.title(
    "CDS Option Value versus Strike",
)

plt.grid(
    True,
)

plt.legend()

plt.show()


# ============================================================================
# 6. CDS OPTION IMPLIED VOLATILITY
# ============================================================================
#
# Each option value above was generated using a constant volatility of 30%.
#
# The implied-volatility routine then solves for the volatility that reproduces
# the calculated option value.
#
# If the pricing and inversion routines are numerically consistent, the
# recovered implied volatility should be close to 30% across the strike range.
#
# This therefore provides a useful numerical consistency test of the CDS
# option implementation.
# ============================================================================

print("\n" + "=" * 78)
print("6. CDS OPTION IMPLIED VOLATILITY")
print("=" * 78)

long_implied_vols = np.asarray(
    long_implied_vols,
)

short_implied_vols = np.asarray(
    short_implied_vols,
)

plt.figure(
    figsize=(9, 6),
)

plt.plot(
    strikes,
    long_implied_vols * 100.0,
    marker="o",
    label="Long Protection",
)

plt.plot(
    strikes,
    short_implied_vols * 100.0,
    marker="o",
    label="Short Protection",
)

plt.axhline(
    volatility * 100.0,
    linestyle="--",
    label="Input Volatility",
)

plt.xlabel(
    "Strike Spread (bp)",
)

plt.ylabel(
    "Implied Volatility (%)",
)

plt.title(
    "CDS Option Implied Volatility versus Strike",
)

plt.grid(
    True,
)

plt.legend()
plt.show()
