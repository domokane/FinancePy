# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 07:54:23 2026

@author: Dominic
"""

import numpy as np

from financepy.utils.date import Date
from financepy.utils.global_vars import ONE_MILLION, CLEAN

from financepy.market.curves.discount_curve import DiscountCurve
from financepy.market.curves.interpolator import InterpTypes
from financepy.market.curves.cds_curve import CDSCurve

from financepy.products.credit.cds import CDS
from financepy.products.credit.cds_fixed_recovery import CDSFixedRecovery


###############################################################################
# Market and contract settings
###############################################################################

value_dt = Date(20, 6, 2026)
step_in_dt = value_dt.add_days(1)
maturity_dt = step_in_dt.next_cds_date(60)

notional = ONE_MILLION
market_spread = 0.0100       # 100 bp
running_cpn = 0.0050         # 50 bp
fixed_recovery = 0.20        # Contractually fixed at 20%

# Flat 3% discount curve
times = np.linspace(0.0, 10.0, 11)
dates = value_dt.add_years(times)
discount_factors = np.exp(-0.03 * times)

discount_curve = DiscountCurve(
    value_dt,
    dates,
    discount_factors,
    InterpTypes.FLAT_FWD_RATES,
)

###############################################################################
# Two contracts with otherwise identical terms
###############################################################################

ordinary_cds = CDS(
    step_in_dt,
    maturity_dt,
    running_cpn,
    notional,
)

fixed_cds = CDSFixedRecovery(
    step_in_dt,
    maturity_dt,
    running_cpn,
    fixed_recovery,
    notional,
)

###############################################################################
# Compare under different issuer curve recovery assumptions
###############################################################################

print("\nCOMPARISON OF ORDINARY AND FIXED-RECOVERY CDS")
print("Market CDS spread: 100 bp")
print("Contract running coupon: 50 bp")
print("Fixed contractual recovery: 20%")
print()

print(
    f"{'Curve Rec':>10} "
    f"{'CDS Par(bp)':>12} "
    f"{'Fixed Par(bp)':>14} "
    f"{'CDS PV':>14} "
    f"{'Fixed PV':>14}"
)

realised_recovery_rates = np.linspace(0.0, 0.6, 61)

for curve_recovery in realised_recovery_rates:

    # The same 100 bp market CDS quote is used for calibration.
    market_cds = CDS(
        step_in_dt,
        maturity_dt,
        market_spread,
    )

    issuer_curve = CDSCurve(
        value_dt,
        [market_cds],
        discount_curve,
        curve_recovery,
    )

    # Ordinary CDS: valuation recovery equals the curve assumption.
    ordinary_par = ordinary_cds.par_spread(
        value_dt,
        issuer_curve,
        curve_recovery,
    )

    ordinary_pv = ordinary_cds.value(
        value_dt,
        issuer_curve,
        curve_recovery,
    )[CLEAN]

    # Fixed-recovery CDS: contractual recovery remains 20%.
    fixed_par = fixed_cds.par_spread(
        value_dt,
        issuer_curve,
        fixed_recovery,
    )

    fixed_pv = fixed_cds.value(
        value_dt,
        issuer_curve,
    )[CLEAN]

    print(
        f"{curve_recovery:10.0%} "
        f"{ordinary_par * 10000:12.2f} "
        f"{fixed_par * 10000:14.2f} "
        f"{ordinary_pv:14,.2f} "
        f"{fixed_pv:14,.2f}"
    )

###############################################################################
# Compare protection payments following default
###############################################################################

print("\nPROTECTION PAYMENT FOLLOWING DEFAULT")
print()

print(
    f"{'Realised Recovery':>20} "
    f"{'Ordinary CDS':>16} "
    f"{'Fixed CDS':>16}"
)


for realised_recovery in [0.0, 0.2, 0.4, 0.6]:

    ordinary_payment = notional * (1.0 - realised_recovery)
    fixed_payment = notional * (1.0 - fixed_recovery)

    print(
        f"{realised_recovery * 100:.0f}\\% & "
        f"{ordinary_payment:,.2f} & "
        f"{fixed_payment:,.2f} \\\\")

###############################################################################
# Recovery sensitivity: curve recovery is bumped, not fixed recovery
###############################################################################

curve_recovery = 0.40

market_cds = CDS(step_in_dt, maturity_dt, market_spread)

issuer_curve = CDSCurve(
    value_dt,
    [market_cds],
    discount_curve,
    curve_recovery,
)

ordinary_recovery01 = ordinary_cds.recovery_dv01(
    value_dt,
    issuer_curve,
    curve_recovery,
)

fixed_recovery01 = fixed_cds.recovery_dv01(
    value_dt,
    issuer_curve,
)

print("\nRECOVERY DV01")
print(f"Ordinary CDS:       {ordinary_recovery01:,.2f}")
print(f"Fixed-recovery CDS: {fixed_recovery01:,.2f}")
