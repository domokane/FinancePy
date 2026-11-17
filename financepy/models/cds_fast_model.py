##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit

STANDARD_RECOVERY_RATE = 0.40
ONE_BP = 0.0001
ONE_PCT = 0.01
KAPPA = 365.0 / 360.0

##############################################################################


@njit(cache=True)
def risky_annuity_flat(t_eff, t_mat, rate, hazard_rate):

    horizon = t_mat - t_eff
    decay = rate + hazard_rate

    x = decay * horizon
    small = np.abs(x) < 1.0e-8
    safe_x = np.where(small, 1.0, x)

    annuity = horizon * np.where(
        small,
        1.0 - x / 2.0 + x * x / 6.0,
        -np.expm1(-x) / safe_x,
    )

    return np.exp(-decay * t_eff) * annuity

##############################################################################


@njit(cache=True)
def prot_leg_pv_flat(
    t_eff,
    t_mat,
    rate,
    hazard_rate,
    recovery_rate,
):

    annuity = risky_annuity_flat(
        t_eff, t_mat, rate, hazard_rate
    )

    return hazard_rate * (1.0 - recovery_rate) * annuity

##############################################################################


@njit(cache=True)
def rpv01_flat(
    t_eff,
    t_mat,
    rate,
    hazard_rate,
    delta,
):

    annuity = risky_annuity_flat(
        t_eff, t_mat, rate, hazard_rate
    )

    clean_rpv01 = KAPPA * annuity
    dirty_rpv01 = clean_rpv01 + delta

    return dirty_rpv01, clean_rpv01

##############################################################################


@njit(cache=True)
def cds_fast_value(
    running_cpn: float,
    spread: float,
    t_eff: float,
    t_mat: float,
    delta: float,
    rate: float,
    curve_recovery: float = STANDARD_RECOVERY_RATE,
    contract_recovery: float = STANDARD_RECOVERY_RATE,
    notional: float = 1000000.0,
    direction: int = 1,
    bump_size: float = ONE_BP,
    rec_bump: float = ONE_PCT,
    fixed_recovery: bool = False
):

    dirty_pv, clean_pv, rpv01_dirty, rpv01_clean = dirty_and_clean_pv(
        running_cpn, spread, rate, t_eff, t_mat, delta,
        curve_recovery, contract_recovery, notional, direction
    )

    spread_up = spread + bump_size
    spread_dn = np.maximum(spread - bump_size, 0.0)

    dirty_pv_spread_up, _, _, _ = dirty_and_clean_pv(
        running_cpn, spread_up, rate, t_eff, t_mat, delta,
        curve_recovery, contract_recovery, notional, direction
    )

    dirty_pv_spread_dn, _, _, _ = dirty_and_clean_pv(
        running_cpn, spread_dn, rate, t_eff, t_mat, delta,
        curve_recovery, contract_recovery, notional, direction
    )

    rate_up = rate + bump_size
    rate_dn = np.maximum(rate - bump_size, 0.0)

    dirty_pv_rate_up, _, _, _ = dirty_and_clean_pv(
        running_cpn, spread, rate_up, t_eff, t_mat, delta,
        curve_recovery, contract_recovery, notional, direction
    )

    dirty_pv_rate_dn, _, _, _ = dirty_and_clean_pv(
        running_cpn, spread, rate_dn, t_eff, t_mat, delta,
        curve_recovery, contract_recovery, notional, direction
    )

    if fixed_recovery:

        dirty_pv_recovery_up, _, _, _ = dirty_and_clean_pv(
            running_cpn, spread, rate, t_eff, t_mat, delta,
            curve_recovery + rec_bump, contract_recovery,
            notional, direction
        )

        dirty_pv_recovery_dn, _, _, _ = dirty_and_clean_pv(
            running_cpn, spread, rate, t_eff, t_mat, delta,
            curve_recovery - rec_bump, contract_recovery,
            notional, direction
        )

    else:

        dirty_pv_recovery_up, _, _, _ = dirty_and_clean_pv(
            running_cpn, spread, rate, t_eff, t_mat, delta,
            curve_recovery + rec_bump, contract_recovery + rec_bump,
            notional, direction
        )

        dirty_pv_recovery_dn, _, _, _ = dirty_and_clean_pv(
            running_cpn, spread, rate, t_eff, t_mat, delta,
            curve_recovery - rec_bump, contract_recovery - rec_bump,
            notional, direction
        )

    spread_dv01 = (
        (dirty_pv_spread_up - dirty_pv_spread_dn)
        * ONE_BP / (spread_up - spread_dn)
    )

    ir_dv01 = (
        (dirty_pv_rate_up - dirty_pv_rate_dn)
        * ONE_BP / (rate_up - rate_dn)
    )

    rec_dv01 = (
        (dirty_pv_recovery_up - dirty_pv_recovery_dn)
        * ONE_PCT / (2.0 * rec_bump)
    )

    return (
        dirty_pv,
        clean_pv,
        rpv01_dirty,
        rpv01_clean,
        spread_dv01,
        ir_dv01,
        rec_dv01,
    )

##############################################################################


@njit(cache=True)
def dirty_and_clean_pv(
    running_cpn: float,
    spread: float,
    rate: float,
    t_eff: float,
    t_mat: float,
    delta: float,
    curve_recovery: float,
    contract_recovery: float,
    notional: float,
    direction: int,
):
    """Dirty and clean CDS PV under flat hazard and discount curves.

    All numerical inputs support NumPy broadcasting.
    """

    horizon = t_mat - t_eff

    hazard_rate = KAPPA * spread / (1.0 - curve_recovery)
    decay_rate = rate + hazard_rate

    x = decay_rate * horizon

    small = np.abs(x) < 1.0e-8

    safe_decay = np.where(small, 1.0, decay_rate)

    annuity_per_df = np.where(
        small,
        horizon * (1.0 - x / 2.0 + x * x / 6.0),
        -np.expm1(-x) / safe_decay,
    )

    risky_annuity = np.exp(-decay_rate * t_eff) * annuity_per_df

    rpv01_clean = KAPPA * risky_annuity
    rpv01_dirty = rpv01_clean + delta

    protection_pv = (
        hazard_rate
        * (1.0 - contract_recovery)
        * risky_annuity
        * notional
    )

    premium_pv_dirty = running_cpn * rpv01_dirty * notional
    premium_pv_clean = running_cpn * rpv01_clean * notional

    dirty_pv = direction * (protection_pv - premium_pv_dirty)
    clean_pv = direction * (protection_pv - premium_pv_clean)

    return dirty_pv, clean_pv, rpv01_dirty, rpv01_clean
