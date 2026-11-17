##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

"""
O'Kane-Turnbull (OT) CDS Valuation Model
--------------------------------------

This module implements a reduced-form valuation model for credit default
swaps based on the methodology described in:

    O'Kane, D. and Turnbull, S. (2003).
    "Valuation of Credit Default Swaps."
    Lehman Brothers, Fixed Income Quantitative Credit Research, April.

The name OT is used to distinguish this implementation from the ISDA
Standard CDS Model. It follows the general valuation methodology of
O'Kane and Turnbull, with subsequent extensions to accommodate modern
CDS market conventions and improve numerical accuracy.

Model assumptions:

    1. Interest rates and default intensities are deterministic.

    2. Discount factors and survival probabilities are interpolated
       using piecewise-constant forward rates and hazard rates.

    3. The recovery rate is constant and independent of default time.

    4. Counterparty credit risk is ignored.

    5. The protection payment is assumed to occur immediately upon
       default.

    6. The protection leg is valued by numerical integration over
       a regular time grid, using analytical integration under
       locally constant hazard and forward interest rates.

    7. The premium leg consists of discrete scheduled coupon payments
       conditional on survival, together with accrued premium payable
       upon default.

    8. Accrued premium on default is approximated in the first coupon
       period and integrated analytically under locally constant
       hazard and forward rates in subsequent coupon periods.

    9. Protection is modelled using the beginning-of-day
       convention. Scheduled premium survival probabilities are
       evaluated at the inclusive accrual end date.
    
    10. The accrued premium rebate is assumed to be exchanged
        on the cash settlement date. Clean risky PV01 is obtained
        by subtracting the present value of this rebate from
        dirty risky PV01.

The protection leg and risky PV01 are calculated independently.
Both clean and dirty risky PV01 values are returned, allowing the
calculation of clean and dirty CDS values, par spreads, upfront
payments and cash settlement amounts.

The numerical implementation is accelerated using Numba.

Although the OT and ISDA models share the same underlying
reduced-form valuation framework and contractual conventions,
they differ in their numerical integration schemes and
approximations for accrued premium payable upon default.

The OT model is not intended to reproduce the ISDA Standard
Model exactly.
"""

import numpy as np
from numba import njit

from ..utils.error import FinError
from ..market.curves.interpolator import InterpTypes, _uinterpolate
from ..utils.global_vars import CLEAN, DIRTY

GLOB_NUM_STEPS_PER_YEAR = 25

USE_FLAT_HAZARD_RATE_INTEGRAL = True

DAY = 1.0 / 365.0

# @njit(
#     fastmath=True,
#     cache=True,
# )
# def rpv01_ot_numba_old(
#     t_eff,
#     accrual_factor_pcd_to_now,
#     payment_times,
#     year_fracs,
#     np_ibor_times,
#     np_ibor_values,
#     np_surv_times,
#     np_surv_values,
# ):
#     """Fast calculation of the risky PV01 of a CDS using NUMBA.
#     The output is a numpy array of the full and clean risky PV01."""

#     method = InterpTypes.FLAT_FWD_RATES.value

#     cpn_accd_indicator = 1

#     # Method 0 : This is the market standard which assumes that the cpn
#     # accrued is treated as though on average default occurs roughly midway
#     # through a cpn period.

#     t_ncd = payment_times[0]

#     # The first cpn is a special case which needs to be handled carefully
#     # taking into account what cpn has already accrued and what has not
#     q_eff = _uinterpolate(t_eff, np_surv_times, np_surv_values, method)
#     q1 = _uinterpolate(t_ncd, np_surv_times, np_surv_values, method)
#     z1 = _uinterpolate(t_ncd, np_ibor_times, np_ibor_values, method)

#     # The accrual factor of the current cpn period, from the previous cpn
#     # date to the next cpn date. The entries of year_fracs are aligned with
#     # payment_times so index 0 is the current period.
#     accrual_factor_ncd = year_fracs[0]

#     # Probability of default between effective date and first coupon payment
#     dq = q_eff - q1

#     dirty_rpv01 = q1 * z1 * accrual_factor_ncd

#     dirty_rpv01 += (
#         z1 * dq * accrual_factor_pcd_to_now
#         * cpn_accd_indicator
#     )

#     dirty_rpv01 += (
#         0.5 * z1 * dq
#         * (accrual_factor_ncd - accrual_factor_pcd_to_now)
#         * cpn_accd_indicator
#     )

#     accrued_adjustment = accrual_factor_pcd_to_now

#     first_accrual_factor = accrual_factor_ncd

#     log_table = []

#     log_table.append([t_ncd, q1, z1, first_accrual_factor, dirty_rpv01])

#     for it in range(1, len(payment_times)):

#         t2 = payment_times[it]
#         q2 = _uinterpolate(t2, np_surv_times, np_surv_values, method)
#         z2 = _uinterpolate(t2, np_ibor_times, np_ibor_values, method)
#         accrual_factor = year_fracs[it]

#         # full cpn is paid at the end of the current period if survives
#         dirty_rpv01 += q2 * z2 * accrual_factor

#         if cpn_accd_indicator == 1:

#             if USE_FLAT_HAZARD_RATE_INTEGRAL:
#                 # TODO This needs to be updated to handle small h+r
#                 tau = accrual_factor
#                 h12 = -np.log(q2 / q1) / tau
#                 r12 = -np.log(z2 / z1) / tau
#                 alpha = h12 + r12
#                 exp_term = 1.0 - np.exp(-alpha * tau) - alpha * tau * np.exp(-alpha * tau)
#                 d_dirty_rpv01 = q1 * z1 * h12 * exp_term / abs(alpha * alpha + 1e-20)
#             else:
#                 d_dirty_rpv01 = 0.50 * (q1 - q2) * z2 * accrual_factor

#             dirty_rpv01 += d_dirty_rpv01

#         log_table.append([t2, q2, z2, accrual_factor, dirty_rpv01])

#         q1 = q2
#         z1 = z2

#     clean_rpv01 = dirty_rpv01 - accrued_adjustment

#     v = np.array([0.0, 0.0])
#     v[DIRTY] = dirty_rpv01
#     v[CLEAN] = clean_rpv01

#     log_table = np.array(log_table)

#     return v, log_table


##############################################################################

@njit(
    fastmath=True,
    cache=True,
)
def rpv01_ot_numba_OLD(
        t_eff,
        accrual_factor_pcd_to_now,
        payment_times,
        year_fracs,
        np_ibor_times,
        np_ibor_values,
        np_surv_times,
        np_surv_values,
        accrual_start_times,
        accrual_end_times,
        settlement_df,
        pay_accrued_on_default=True):
    """Calculate dirty and clean risky PV01 using the OT methodology.

    Scheduled premium payments are discounted to their payment dates.
    Accrual boundaries determine the default-accrual integration periods.

    Accrued premium on default is approximated in the first coupon
    period and integrated analytically under locally constant hazard
    and forward interest rates in subsequent coupon periods.

    All times are measured on an ACT/365F basis. The contractual
    accrual fractions are supplied separately in year_fracs.

    Accrual end times represent exclusive accrual boundaries.
    """

    method = InterpTypes.FLAT_FWD_RATES.value

    cpn_accd_indicator = 1

    num_payments = len(payment_times)

    dirty_rpv01 = 0.0

    log_table = np.zeros((num_payments, 5))

    for it in range(num_payments):

        t_pay = payment_times[it]
        t_start = accrual_start_times[it]
        t_end = accrual_end_times[it]

        accrual_factor = year_fracs[it]

        # Survival is measured at the end of the covered accrual
        # period, whereas discounting uses the payment date.
        q_end = _uinterpolate(
            t_end,
            np_surv_times,
            np_surv_values,
            method,
        )

        z_pay = _uinterpolate(
            t_pay,
            np_ibor_times,
            np_ibor_values,
            method,
        )

        # Scheduled coupon conditional on survival.
        dirty_rpv01 += q_end * z_pay * accrual_factor

        t_obs = t_end - DAY

        accrual_factor = year_fracs[it]

        q_end = _uinterpolate(
            t_obs,
            np_surv_times,
            np_surv_values,
            method,
        )

        if pay_accrued_on_default:

            # Default accrual cannot start before protection begins.

            t0 = max(0.0, t_start - DAY, t_eff)

            tau = t_obs - t0

            if tau > 0.0:

                q0 = _uinterpolate(
                    t0,
                    np_surv_times,
                    np_surv_values,
                    method,
                )

                dq = q0 - q_end

                if it == 0:

                    # First coupon period: retain the OT midpoint
                    # approximation, including premium already accrued.

                    accrued = accrual_factor_pcd_to_now

                    remaining = accrual_factor - accrued

                    d_dirty_rpv01 = (
                        z_pay * dq
                        * (accrued + 0.5 * remaining)
                    )

                elif USE_FLAT_HAZARD_RATE_INTEGRAL:

                    # Subsequent periods: analytical integration
                    # assuming constant hazard and forward rates
                    # over the accrual interval.

                    z0 = _uinterpolate(
                        t0,
                        np_ibor_times,
                        np_ibor_values,
                        method,
                    )

                    z_end = _uinterpolate(
                        t_end,
                        np_ibor_times,
                        np_ibor_values,
                        method,
                    )

                    h12 = -np.log(q_end / q0) / tau
                    r12 = -np.log(z_end / z0) / tau

                    alpha = h12 + r12
                    x = alpha * tau

                    # Convert calendar time to contractual
                    # premium accrual time.
                    kappa = accrual_factor / (t_end - t_start)

                    # First exponential moment:
                    # integral_0^tau u * exp(-alpha*u) du

                    if abs(x) < 1.0e-4:

                        integral = tau * tau * (
                            0.5
                            - x / 3.0
                            + x * x / 8.0
                            - x * x * x / 30.0
                        )

                    else:

                        integral = (
                            -np.expm1(-x)
                            - x * np.exp(-x)
                        ) / (alpha * alpha)

                    d_dirty_rpv01 = (
                        kappa * q0 * z0 * h12 * integral
                    )

                else:

                    # Simple midpoint approximation.
                    d_dirty_rpv01 = (
                        0.5 * dq * z_pay * accrual_factor
                    )

                dirty_rpv01 += d_dirty_rpv01

        log_table[it, 0] = t_pay
        log_table[it, 1] = q_end
        log_table[it, 2] = z_pay
        log_table[it, 3] = accrual_factor
        log_table[it, 4] = dirty_rpv01

    # Accrued premium rebate is paid on the cash settlement date.
    clean_rpv01 = dirty_rpv01 - settlement_df * accrual_factor_pcd_to_now

    v = np.zeros(2)

    v[DIRTY] = dirty_rpv01
    v[CLEAN] = clean_rpv01

    return v, log_table

########################################################################################


@njit(fastmath=True, cache=True)
def rpv01_ot_numba(
    t_eff,
    accrual_factor_pcd_to_now,
    payment_times,
    year_fracs,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    accrual_start_times,
    accrual_end_times,
    settlement_df,
    pay_accrued_on_default=True,
):
    """Calculate dirty and clean risky PV01 using the OT methodology.

    Scheduled premium payments are discounted to their payment dates.
    Survival is measured at the inclusive accrual end date.

    Accrued premium on default is approximated in the first coupon
    period and integrated analytically under locally constant hazard
    and forward rates in subsequent coupon periods.

    All times are measured on an ACT/365F basis. Contractual accrual
    fractions are supplied separately in year_fracs.

    Accrual end times represent exclusive accrual boundaries.

    The effective time is the protection integration start time,
    already adjusted for beginning-of-day protection.

    The accrued premium rebate is discounted from the cash
    settlement date to the valuation date.
    """

    method = InterpTypes.FLAT_FWD_RATES.value

    num_payments = len(payment_times)

    dirty_rpv01 = 0.0

    log_table = np.zeros((num_payments, 5))

    for it in range(num_payments):

        t_pay = payment_times[it]
        t_start = accrual_start_times[it]
        t_end = accrual_end_times[it]

        accrual_factor = year_fracs[it]

        # Beginning-of-day protection convention.
        t_obs = t_end - DAY
        t_accrual_start = t_start - DAY

        # Survival is measured at the observation date.
        q_end = _uinterpolate(
            t_obs,
            np_surv_times,
            np_surv_values,
            method,
        )

        # Discounting uses the actual payment date.
        z_pay = _uinterpolate(
            t_pay,
            np_ibor_times,
            np_ibor_values,
            method,
        )

        # Scheduled coupon conditional on survival.
        dirty_rpv01 += q_end * z_pay * accrual_factor

        if pay_accrued_on_default:

            # Default accrual cannot start before protection.
            t0 = max(0.0, t_accrual_start, t_eff)

            tau = t_obs - t0

            if tau > 0.0:

                q0 = _uinterpolate(
                    t0,
                    np_surv_times,
                    np_surv_values,
                    method,
                )

                dq = q0 - q_end

                if it == 0:

                    # OT midpoint approximation for the first period,
                    # including premium already accrued.

                    accrued = accrual_factor_pcd_to_now

                    remaining = accrual_factor - accrued

                    d_dirty_rpv01 = (
                        z_pay
                        * dq
                        * (accrued + 0.5 * remaining)
                    )

                elif USE_FLAT_HAZARD_RATE_INTEGRAL:

                    # Analytical integration under locally constant
                    # hazard and forward interest rates.

                    z0 = _uinterpolate(
                        t0,
                        np_ibor_times,
                        np_ibor_values,
                        method,
                    )

                    z_end = _uinterpolate(
                        t_obs,
                        np_ibor_times,
                        np_ibor_values,
                        method,
                    )

                    h12 = -np.log(q_end / q0) / tau
                    r12 = -np.log(z_end / z0) / tau

                    alpha = h12 + r12

                    x = alpha * tau

                    # Convert calendar time to contractual
                    # premium accrual time.
                    kappa = accrual_factor / (t_end - t_start)

                    # First exponential moment:
                    # integral_0^tau u * exp(-alpha*u) du

                    if abs(x) < 1.0e-4:

                        integral = tau * tau * (
                            0.5
                            - x / 3.0
                            + x * x / 8.0
                            - x * x * x / 30.0
                        )

                    else:

                        integral = (
                            -np.expm1(-x)
                            - x * np.exp(-x)
                        ) / (alpha * alpha)

                    d_dirty_rpv01 = (
                        kappa
                        * q0
                        * z0
                        * h12
                        * integral
                    )

                else:

                    # Simple midpoint approximation.
                    d_dirty_rpv01 = (
                        0.5 * dq * z_pay * accrual_factor
                    )

                dirty_rpv01 += d_dirty_rpv01

        log_table[it, 0] = t_pay
        log_table[it, 1] = q_end
        log_table[it, 2] = z_pay
        log_table[it, 3] = accrual_factor
        log_table[it, 4] = dirty_rpv01

    # Accrued premium rebate is paid on the cash settlement date.
    clean_rpv01 = (
        dirty_rpv01
        - settlement_df * accrual_factor_pcd_to_now
    )

    v = np.zeros(2)

    v[DIRTY] = dirty_rpv01
    v[CLEAN] = clean_rpv01

    return v, log_table

########################################################################################


@njit(
    fastmath=True,
    cache=True,
)
def prot_leg_pv_ot_numba(
    t_eff,
    t_mat,
    np_ibor_times,
    np_ibor_values,
    np_surv_times,
    np_surv_values,
    recovery_rate,
):
    """Fast calculation of the CDS protection leg PV using NUMBA to speed up
    the numerical integration over time."""

    num_steps_per_year = GLOB_NUM_STEPS_PER_YEAR

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

    if USE_FLAT_HAZARD_RATE_INTEGRAL:

        log_z1 = np.log(z1)
        log_q1 = np.log(q1)

        for _ in range(0, num_steps):
            t = t + dt
            z2 = _uinterpolate(t, np_ibor_times, np_ibor_values, method)
            q2 = _uinterpolate(t, np_surv_times, np_surv_values, method)

            log_z2 = np.log(z2)
            log_q2 = np.log(q2)

            h12 = -(log_q2 - log_q1) / dt
            r12 = -(log_z2 - log_z1) / dt

            alpha = h12 + r12
            x = alpha * dt

            if abs(x) < 1e-8:
                integral = dt * (
                    1.0 - x / 2.0 + x * x / 6.0
                )
            else:
                integral = -np.expm1(-x) / alpha

            dprot_pv = h12 * q1 * z1 * integral

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
