##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

""" ISDA-convention CDS leg valuation (v1.8.3 conventions).

Independent implementation of the reduced-form leg integrals, with the
original half-day accrual-on-default convention used by ISDA CDS Standard
Model v1.8.3. This module is not the official ISDA software and does not
implement the standard-model curve builders or spread/upfront converter.

DISCLAIMER:
This CDS Model is not guaranteed to reproduce the ISDA Standard Model exactly.

Assumptions
-----------
* ACT/365F curve times, relative to the common valuation-date anchor.
* Log-linear discount factors and survival probabilities; flat forward
  interest rates and hazard rates between their respective nodes.
* Standard start-of-day protection, including the maturity day. The
  pricing step-in date is at least one calendar day after the valuation
  date (end-of-day valuation).
* Constant recovery, immediate protection payment at default, no
  counterparty risk.
* Discrete scheduled premiums discounted to their payment dates, with
  survival measured to the day before the exclusive accrual end date.
* Default accrual integrated over the combined curve grid with the
  original ISDA half-day bias. ACT/360 is the standard premium day count;
  the supplied year fractions must be consistent with the accrual dates.
* Clean RPV01 removes the step-in accrued fraction discounted to the
  cash settlement date. All outputs are valuation-date PVs, not settlement
  prices. Divide by the settlement discount factor to obtain the latter.

The caller supplies dates, accrual fractions and discount/survival curves.
These kernels have no dependencies on CDS, CDSCurve or calendar classes.
They return protection PV per unit notional (including loss given default)
and a dirty/clean RPV01 pair, with a five-column coupon valuation table.

Reference conventions: https://www.cdsmodel.com/ (version 1.8.3).
This is a compatibility candidate, not a claim of official certification.
"""

import numpy as np
from numba import njit

DAY = 1.0 / 365.0
HALF_DAY = 0.5 / 365.0
TIME_TOL = 1.0e-12


@njit(cache=True)
def _curve_logs(times, factors):
    """Validate a one-dimensional curve and return logarithmic factors."""
    n = len(times)
    if n == 0 or n != len(factors):
        raise ValueError("Curve times and factors must have equal nonzero length.")
    if times[0] < 0.0 or times[-1] <= 0.0:
        raise ValueError("Curve requires at least one positive-time node.")
    logs = np.empty(n, dtype=np.float64)
    for i in range(n):
        if not np.isfinite(times[i]) or not np.isfinite(factors[i]):
            raise ValueError("Curve times and factors must be finite.")
        if factors[i] <= 0.0:
            raise ValueError("Discount and survival factors must be positive.")
        if i > 0 and times[i] <= times[i - 1]:
            raise ValueError("Curve times must be strictly increasing.")
        logs[i] = np.log(factors[i])
    if times[0] == 0.0 and abs(logs[0]) > TIME_TOL:
        raise ValueError("A curve anchored at time zero must have factor one.")
    return logs


@njit(cache=True)
def _log_factor(t, times, logs):
    """Log-linear interpolation with flat-forward extrapolation."""
    if t == 0.0:
        return 0.0
    n = len(times)
    if t <= times[0]:
        j = 1 if times[0] == 0.0 else 0
        return logs[j] * t / times[j]
    j = np.searchsorted(times, t)
    if j == 0:
        return logs[0] * t / times[0]
    if j == n:
        if n == 1:
            return logs[0] * t / times[0]
        j = n - 1
    i = j - 1
    w = (t - times[i]) / (times[j] - times[i])
    return logs[i] + w * (logs[j] - logs[i])


@njit(cache=True)
def _exp_moments(x):
    """Integrals of exp(-x*u) and u*exp(-x*u) over u in [0, 1]."""
    if abs(x) < 0.01:
        # Taylor coefficients (-1)^n / (n! * (n + k + 1)), k=0,1.
        p0 = 1.0 + x * (-0.5 + x * (1.0 / 6.0 + x * (
            -1.0 / 24.0 + x * (1.0 / 120.0 + x * (
                -1.0 / 720.0 + x / 5040.0)))))
        p1 = 0.5 + x * (-1.0 / 3.0 + x * (1.0 / 8.0 + x * (
            -1.0 / 30.0 + x * (1.0 / 144.0 + x * (
                -1.0 / 840.0 + x / 5760.0)))))
    else:
        em1 = np.expm1(-x)
        p0 = -em1 / x
        p1 = (-em1 - x * np.exp(-x)) / (x * x)
    return p0, p1


@njit(cache=True)
def _default_integrals(start, end, accrual_origin, nodes,
                       discount_times, discount_logs,
                       survival_times, survival_logs):
    """Unit default PV and time-weighted default PV on [start, end].

    accrual_origin is the zero of the time weight, including any bias.
    The combined grid has already been prepared by the caller.
    """
    if end <= start + TIME_TOL:
        return 0.0, 0.0

    u = start
    log_z0 = _log_factor(u, discount_times, discount_logs)
    log_q0 = _log_factor(u, survival_times, survival_logs)
    j = np.searchsorted(nodes, start + TIME_TOL)
    default_pv = 0.0
    accrual_pv = 0.0

    while u < end:
        if j < len(nodes) and nodes[j] < end - TIME_TOL:
            v = nodes[j]
            j += 1
        else:
            v = end

        dt = v - u
        log_z1 = _log_factor(v, discount_times, discount_logs)
        log_q1 = _log_factor(v, survival_times, survival_logs)

        # Integrated hazard and forward interest rate over the interval.
        dh = log_q0 - log_q1
        df = log_z0 - log_z1
        p0, p1 = _exp_moments(dh + df)

        weight = dh * np.exp(log_z0 + log_q0)
        default_pv += weight * p0
        accrual_pv += weight * ((u - accrual_origin) * p0 + dt * p1)

        u = v
        log_z0 = log_z1
        log_q0 = log_q1

    return default_pv, accrual_pv


@njit(cache=True)
def prot_leg_pv_isda_numba(t_stepin, t_mat, discount_times, discount_values,
                           survival_times, survival_values, recovery_rate):
    """Protection PV per unit notional, including loss given default.

    t_stepin: pricing step-in date, in ACT/365F years from valuation date.
              Pass the actual (unshifted) date, not step-in minus one day.
    t_mat: unadjusted contractual maturity date, in ACT/365F years.

    Start-of-day protection is implemented by integrating from step-in
    minus one calendar day to the end of the maturity day.
    """
    if not np.isfinite(t_stepin) or not np.isfinite(t_mat):
        raise ValueError("Protection dates must be finite.")
    if t_stepin < DAY - TIME_TOL:
        raise ValueError("Pricing step-in must be at least valuation date + 1 day.")
    if not 0.0 <= recovery_rate <= 1.0:
        raise ValueError("Contract recovery must lie in [0, 1].")
    if t_stepin > t_mat + TIME_TOL:
        return 0.0

    discount_logs = _curve_logs(discount_times, discount_values)
    survival_logs = _curve_logs(survival_times, survival_values)
    nodes = np.unique(np.concatenate((discount_times, survival_times)))

    start = max(0.0, t_stepin - DAY)
    default_pv, _ = _default_integrals(
        start, t_mat, start, nodes,
        discount_times, discount_logs, survival_times, survival_logs,
    )
    return (1.0 - recovery_rate) * default_pv


@njit(cache=True)
def rpv01_isda_numba(t_stepin, accrued_fraction, payment_times, year_fracs,
                     discount_times, discount_values, survival_times,
                     survival_values, accrual_start_times, accrual_end_times,
                     settlement_df, pay_accrued_on_default=True):
    """Return (RPV01 array, coupon log table), both as valuation-date PVs.

    The RPV01 array contains [dirty, clean], consistent with FinancePy.
    t_stepin is the actual pricing step-in date (unshifted).
    payment_times are actual business-day-adjusted payment dates.
    accrual_start_times are inclusive contractual accrual start dates.
    accrual_end_times are EXCLUSIVE accrual end dates. For the last coupon,
    this is the unadjusted maturity date PLUS one calendar day. For earlier
    coupons it is the next accrual start date. Do not substitute payment
    times for these dates, and do not shift them back a day before calling.
    year_fracs are contractual day-count fractions for the full periods.
    accrued_fraction is the day-count fraction accrued to pricing step-in.
    settlement_df is the discount factor from today to cash settlement.

    Table columns: payment time, survival at observation date, payment DF,
    full accrual fraction, cumulative dirty RPV01. Expired rows contribute
    zero and are omitted from the table. Normally the CDS wrapper already
    removes expired coupon periods before calling.
    """
    if not np.isfinite(t_stepin) or t_stepin < DAY - TIME_TOL:
        raise ValueError("Pricing step-in must be at least valuation date + 1 day.")
    if not np.isfinite(accrued_fraction) or accrued_fraction < 0.0:
        raise ValueError("Accrued fraction must be finite and non-negative.")
    if not np.isfinite(settlement_df) or settlement_df <= 0.0:
        raise ValueError("Settlement discount factor must be positive and finite.")

    n = len(payment_times)
    if (len(year_fracs) != n or len(accrual_start_times) != n
            or len(accrual_end_times) != n):
        raise ValueError("Payment times, accrual boundaries and fractions must align.")
    table = np.empty((n, 5), dtype=np.float64)
    if n == 0 or t_stepin > accrual_end_times[-1] - DAY + TIME_TOL:
        return np.zeros(2, dtype=np.float64), table[:0]

    discount_logs = _curve_logs(discount_times, discount_values)
    survival_logs = _curve_logs(survival_times, survival_values)
    nodes = np.unique(np.concatenate((discount_times, survival_times)))
    protection_start = max(0.0, t_stepin - DAY)
    dirty = 0.0
    k = 0

    for i in range(n):
        a = accrual_start_times[i]
        b = accrual_end_times[i]
        p = payment_times[i]
        fraction = year_fracs[i]
        if (not np.isfinite(a) or not np.isfinite(b) or not np.isfinite(p)
                or not np.isfinite(fraction) or b <= a or fraction < 0.0):
            raise ValueError("Invalid coupon accrual period.")
        if b <= t_stepin + TIME_TOL:
            continue

        # Survival is observed at end of the last accrued calendar day.
        a_obs = a - DAY
        b_obs = b - DAY
        q = np.exp(_log_factor(b_obs, survival_times, survival_logs))
        z = np.exp(_log_factor(p, discount_times, discount_logs))
        dirty += fraction * q * z

        if pay_accrued_on_default:
            start = max(protection_start, a_obs)
            _, weighted = _default_integrals(
                start, b_obs, a_obs - HALF_DAY, nodes,
                discount_times, discount_logs, survival_times, survival_logs,
            )
            # Full-period year fraction per unit of ACT/365F time.
            dirty += fraction / (b - a) * weighted

        table[k, 0] = p
        table[k, 1] = q
        table[k, 2] = z
        table[k, 3] = fraction
        table[k, 4] = dirty
        k += 1

    clean = dirty - accrued_fraction * settlement_df
    return np.array([dirty, clean]), table[:k]
