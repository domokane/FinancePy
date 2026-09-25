#
# Copyright (C) 2018-2026 Dominic O'Kane
#

import numpy as np
import os
from financepy.market.curves.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_deposit import IborDeposit
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.global_types import SwapTypes
from financepy.products.credit.cds import CDS
from financepy.market.curves.cds_curve import CDSCurve


def build_ibor_curve(value_dt: Date, bump=0.0):

    m = 1.0

    settle_dt = value_dt
    dc_type = DayCountTypes.ACT_360

    spot_days = 0
    settle_dt = value_dt.add_days(spot_days)

    maturity_dt = settle_dt.add_months(1)
    depo1 = IborDeposit(settle_dt, maturity_dt, m * 0.0016, dc_type)

    maturity_dt = settle_dt.add_months(2)
    depo2 = IborDeposit(settle_dt, maturity_dt, m * 0.0020, dc_type)

    maturity_dt = settle_dt.add_months(3)
    depo3 = IborDeposit(settle_dt, maturity_dt, m * 0.0024, dc_type)

    maturity_dt = settle_dt.add_months(6)
    depo4 = IborDeposit(settle_dt, maturity_dt, m * 0.0033, dc_type)

    maturity_dt = settle_dt.add_months(12)
    depo5 = IborDeposit(settle_dt, maturity_dt, m * 0.0056, dc_type)

    depos = []
    depos.append(depo1)
    depos.append(depo2)
    depos.append(depo3)
    depos.append(depo4)
    depos.append(depo5)

    fras = []

    spot_days = 2
    settle_dt = value_dt.add_days(spot_days)

    swaps = []
    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    fixed_freq = FrequencyTypes.SEMI_ANNUAL

    maturity_dt = settle_dt.add_months(24)
    swap1 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0044 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap1)

    maturity_dt = settle_dt.add_months(36)
    swap2 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0078 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap2)

    maturity_dt = settle_dt.add_months(48)
    swap3 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0119 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap3)

    maturity_dt = settle_dt.add_months(60)
    swap4 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0158 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap4)

    maturity_dt = settle_dt.add_months(72)
    swap5 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0192 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap5)

    maturity_dt = settle_dt.add_months(84)
    swap6 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0219 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap6)

    maturity_dt = settle_dt.add_months(96)
    swap7 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0242 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap7)

    maturity_dt = settle_dt.add_months(108)
    swap8 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0261 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap8)

    maturity_dt = settle_dt.add_months(120)
    swap9 = IborSwap(
        settle_dt,
        maturity_dt,
        SwapTypes.PAY,
        m * 0.0276 + bump,
        fixed_freq,
        dc_type,
    )
    swaps.append(swap9)

    libor_curve = IborSingleCurve(value_dt, depos, fras, swaps)

    return libor_curve


##############################################################################

# ============================================================================
# SUPPORTING FUNCTION - BUILD ISSUER CREDIT CURVE
# ============================================================================
#
# We use a simple CDS term structure to construct the issuer curve.
#
# The optional spread_bump parameter moves every CDS calibration spread by
# the same amount.
#
# This allows us to test spread_dv01() independently.
#
# IMPORTANT:
#
# When recovery changes, the issuer curve must be recalibrated because the
# hazard rates implied by market CDS spreads depend upon recovery.
# ============================================================================


def build_issuer_curve(trade_dt, step_in_dt, libor_curve, recovery_rate, spd_bump=0.0, rec_bump=0.0):
    """Build a simple flat issuer credit curve."""

    value_dt = trade_dt

    cds_mkt_contracts = []
    cds_cpn = 0.005743 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(6)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.007497 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(12)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.011132 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(24)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.013932 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(36)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.015764 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(48)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.017366 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(60)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.020928 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(84)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    cds_cpn = 0.022835 + spd_bump
    maturity_dt = step_in_dt.next_cds_date(120)
    cds = CDS(step_in_dt, maturity_dt, cds_cpn)
    cds_mkt_contracts.append(cds)

    issuer_curve = CDSCurve(
        value_dt,
        cds_mkt_contracts,
        libor_curve,
        recovery_rate + rec_bump,
    )

    return issuer_curve


##############################################################################


def load_heterogeneous_issuer_curves(
    value_dt,
    step_in_dt,
    libor_curve,
):

    maturity_3yr = value_dt.next_cds_date(36)
    maturity_5yr = value_dt.next_cds_date(60)
    maturity_7yr = value_dt.next_cds_date(84)
    maturity_10yr = value_dt.next_cds_date(120)

    path = os.path.join(
        os.path.dirname(__file__),
        "data",
        "CDX_NA_IG_S7_SPREADS.csv",
    )

    issuer_curves = []

    with open(path, "r") as data_file:
        data = data_file.readlines()

    for row in data[1:]:

        split_row = row.split(",")

        spd_3yr = float(split_row[1]) / 10000.0
        spd_5yr = float(split_row[2]) / 10000.0
        spd_7yr = float(split_row[3]) / 10000.0
        spd_10yr = float(split_row[4]) / 10000.0
        recovery_rate = float(split_row[5])

        cds_3yr = CDS(
            step_in_dt,
            maturity_3yr,
            spd_3yr,
        )

        cds_5yr = CDS(
            step_in_dt,
            maturity_5yr,
            spd_5yr,
        )

        cds_7yr = CDS(
            step_in_dt,
            maturity_7yr,
            spd_7yr,
        )

        cds_10yr = CDS(
            step_in_dt,
            maturity_10yr,
            spd_10yr,
        )

        cds_contracts = [
            cds_3yr,
            cds_5yr,
            cds_7yr,
            cds_10yr,
        ]

        issuer_curve = CDSCurve(
            value_dt,
            cds_contracts,
            libor_curve,
            recovery_rate,
        )

        issuer_curves.append(issuer_curve)

    return issuer_curves


# ============================================================================
# SUPPORTING FUNCTION - HOMOGENEOUS CREDIT CURVES
# ============================================================================
#
# Construct identical CDS curves for every name in the basket.
#
# A homogeneous basket is useful because it isolates the effect of:
#
#   - default order
#   - correlation
#   - copula choice
#
# without introducing differences in the individual issuer spread curves.
# ============================================================================


def build_homogeneous_issuer_curves(
    value_dt,
    libor_curve,
    cds_spd_3yr,
    cds_spd_5yr,
    cds_spd_7yr,
    cds_spd_10yr,
    num_credits
):

    step_in_dt = value_dt.add_days(1)
    maturity_3yr = step_in_dt.next_cds_date(36)
    maturity_5yr = step_in_dt.next_cds_date(60)
    maturity_7yr = step_in_dt.next_cds_date(84)
    maturity_10yr = step_in_dt.next_cds_date(120)

    recovery_rate = 0.40

    cds_3yr = CDS(
        value_dt,
        maturity_3yr,
        cds_spd_3yr,
    )

    cds_5yr = CDS(
        value_dt,
        maturity_5yr,
        cds_spd_5yr,
    )

    cds_7yr = CDS(
        value_dt,
        maturity_7yr,
        cds_spd_7yr,
    )

    cds_10yr = CDS(
        value_dt,
        maturity_10yr,
        cds_spd_10yr,
    )

    contracts = [
        cds_3yr,
        cds_5yr,
        cds_7yr,
        cds_10yr,
    ]

    issuer_curve = CDSCurve(
        value_dt,
        contracts,
        libor_curve,
        recovery_rate,
    )

    issuer_curves = []

    for _ in range(num_credits):

        issuer_curves.append(issuer_curve)

    return issuer_curves


###############################################################################

# ============================================================================
# SUPPORTING FUNCTION - CHANGE DISCOUNT CURVE BUT FREEZE CREDIT CURVE
# ============================================================================
#
# This helper is used to understand CDS interest-rate risk.
#
# A CDSCurve contains both:
#
#       1. the interest-rate discount curve
#       2. the calibrated survival probabilities
#
# Normally, after changing interest rates, we recalibrate the CDS curve.
#
# For explanatory purposes we also want to answer:
#
#       What happens if interest rates change but the survival
#       probabilities are held fixed?
#
# We therefore construct an empty CDSCurve using the bumped interest-rate
# curve and then copy the survival-probability term structure from the
# original issuer curve.
#
# This curve is NOT a recalibrated market curve. It is used only to decompose
# the interest-rate sensitivity into:
#
#       direct discounting effect
#
# and
#
#       credit-curve recalibration effect.
#
# ============================================================================


def build_frozen_issuer_curve(
    base_issuer_curve,
    bumped_libor_curve,
    value_dt,
    recovery_rate,
):

    frozen_curve = CDSCurve(
        value_dt,
        [],
        bumped_libor_curve,
        recovery_rate,
    )

    # Copy the calibrated credit term structure.
    #
    # _times contains the credit-curve times.
    # _qs contains the corresponding survival probabilities.

    frozen_curve._times = np.array(
        base_issuer_curve._times,
        copy=True,
    )

    frozen_curve._qs = np.array(
        base_issuer_curve._qs,
        copy=True,
    )

    return frozen_curve
