##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64, int64

from ..utils.error import FinError
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.global_vars import CLEAN, DIRTY

########################################################################################
# TODO: Perform protection leg pv analytically using fact that hazard rate and
#       interest rates are flat between their combined node points. Right now I
#       do not find the protection leg PV calculations to be a bottleneck,
#       especially given the speedup benefits of using NUMBA.
########################################################################################

USE_FLAT_HAZARD_RATE_INTEGRAL = True
STANDARD_RECOVERY_RATE = 0.40
GLOB_NUM_STEPS_PER_YEAR = 25
ONE_BP = 0.0001
ONE_PCT = 0.01

# Premium accrues on an Act/360 basis while time is measured in
# calendar years, hence the 365/360 adjustment to the clean RPV01.
KAPPA = 365.0 / 360.0


@njit(
    fastmath=True,
    cache=True,
)
def rpv01_numba(
    t_eff,
    accrual_factor_pcd_to_now,
    payment_times,
    year_fracs,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    rpv01_method
):
    """Fast calculation of the risky PV01 of a CDS using NUMBA.
    The output is a numpy array of the full and clean risky PV01."""

    method = InterpTypes.FLAT_FWD_RATES.value

    cpn_accd_indicator = 1

    # Method 0 : This is the market standard which assumes that the cpn
    # accrued is treated as though on average default occurs roughly midway
    # through a cpn period.

    t_ncd = payment_times[0]

    # The first cpn is a special case which needs to be handled carefully
    # taking into account what cpn has already accrued and what has not
    q_eff = _uinterpolate(t_eff, np_surv_times, np_surv_values, method)
    q1 = _uinterpolate(t_ncd, np_surv_times, np_surv_values, method)
    z1 = _uinterpolate(t_ncd, np_ibor_times, np_ibor_values, method)

    # The accrual factor of the current cpn period, from the previous cpn
    # date to the next cpn date. The entries of year_fracs are aligned with
    # payment_times so index 0 is the current period.
    accrual_factor_ncd = year_fracs[0]

    # Probability of default between effective date and first coupon payment
    dq = q_eff - q1

    dirty_rpv01 = q1 * z1 * accrual_factor_ncd

    dirty_rpv01 += (
        z1 * dq * accrual_factor_pcd_to_now
        * cpn_accd_indicator
    )

    dirty_rpv01 += (
        0.5 * z1 * dq
        * (accrual_factor_ncd - accrual_factor_pcd_to_now)
        * cpn_accd_indicator
    )

    accrued_adjustment = accrual_factor_pcd_to_now

    first_accrual_factor = accrual_factor_ncd

    log_table = []

    log_table.append([t_ncd, q1, z1, first_accrual_factor, dirty_rpv01])

    for it in range(1, len(payment_times)):

        t2 = payment_times[it]
        q2 = _uinterpolate(t2, np_surv_times, np_surv_values, method)
        z2 = _uinterpolate(t2, np_ibor_times, np_ibor_values, method)
        accrual_factor = year_fracs[it]

        # full cpn is paid at the end of the current period if survives
        dirty_rpv01 += q2 * z2 * accrual_factor

        if cpn_accd_indicator == 1:

            if USE_FLAT_HAZARD_RATE_INTEGRAL:
                # TODO This needs to be updated to handle small h+r
                tau = accrual_factor
                h12 = -np.log(q2 / q1) / tau
                r12 = -np.log(z2 / z1) / tau
                alpha = h12 + r12
                exp_term = 1.0 - np.exp(-alpha * tau) - alpha * tau * np.exp(-alpha * tau)
                d_dirty_rpv01 = q1 * z1 * h12 * exp_term / abs(alpha * alpha + 1e-20)
            else:
                d_dirty_rpv01 = 0.50 * (q1 - q2) * z2 * accrual_factor

            dirty_rpv01 += d_dirty_rpv01

        log_table.append([t2, q2, z2, accrual_factor, dirty_rpv01])

        q1 = q2
        z1 = z2

    clean_rpv01 = dirty_rpv01 - accrued_adjustment

    v = np.array([0.0, 0.0])
    v[DIRTY] = dirty_rpv01
    v[CLEAN] = clean_rpv01

    log_table = np.array(log_table)

    return v, log_table

########################################################################################


@njit(
    fastmath=True,
    cache=True,
)
def prot_leg_pv_numba(
    t_eff,
    t_mat,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    recovery_rate,
    num_steps_per_year,
    prot_method,
):
    """Fast calculation of the CDS protection leg PV using NUMBA to speed up
    the numerical integration over time."""

    method = InterpTypes.FLAT_FWD_RATES.value

    if t_eff < 0.0:
        raise FinError("Error: Protection leg starts in past: t_eff < 0")

    if t_mat <= t_eff:
        return 0.0

    if num_steps_per_year <= 0:
        raise ValueError("num_steps_per_year must be positive")

    num_steps = max(
        1,
        int((t_mat - t_eff) * num_steps_per_year + 0.50),
    )

    dt = (t_mat - t_eff) / num_steps

    t = t_eff
    z1 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
    q1 = _uinterpolate(t, np_surv_times, np_surv_values, method)

    prot_pv = 0.0
    small = 1e-8

    if USE_FLAT_HAZARD_RATE_INTEGRAL:

        log_z1 = np.log(z1)
        log_q1 = np.log(q1)

        for _ in range(0, num_steps):
            t = t + dt
            z2 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
            q2 = _uinterpolate(t, np_surv_times, np_surv_values, method)

            log_z2 = np.log(z2)
            log_q2 = np.log(q2)

            # This needs to be updated to handle small h+r
            # h12 = -log(q2 / q1) / dt
            # r12 = -log(z2 / z1) / dt

            h12 = -(log_q2 - log_q1) / dt
            r12 = -(log_z2 - log_z1) / dt

#            exp_term = np.exp(-(r12 + h12) * dt)
#            dprot_pv = h12 * (1.0 - exp_term) * q1 * z1 / (abs(h12 + r12) + small)

            ##################################################################
            alpha = h12 + r12
            x = alpha * dt

            if abs(x) < 1e-8:
                integral = dt * (
                    1.0 - x / 2.0 + x * x / 6.0
                )
            else:
                integral = -np.expm1(-x) / alpha

            dprot_pv = h12 * q1 * z1 * integral
            ##################################################################

            prot_pv += dprot_pv

            q1, z1 = q2, z2
            log_q1, log_z1 = log_q2, log_z2

    else:

        for _ in range(0, num_steps):
            t += dt
            z2 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
            q2 = _uinterpolate(t, np_surv_times, np_surv_values, method)
            dq = q1 - q2
            dprot_pv = 0.5 * (z1 + z2) * dq
            prot_pv += dprot_pv
            q1 = q2
            z1 = z2

    prot_pv = prot_pv * (1.0 - recovery_rate)
    return prot_pv


########################################################################################
