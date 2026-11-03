##############################################################################
# Copyright (C) 2018-2026 Dominic O'Kane
##############################################################################

import numpy as np

from numba import njit


@njit(fastmath=True, cache=True)
def default_times_from_uniforms_numba(u1,
                                      u2,
                                      curve_times,
                                      curve_qs):

    num_credits, num_trials = u1.shape

    corr_times = np.empty(
        (num_credits, 2 * num_trials)
    )

    for i_credit in range(num_credits):

        times = curve_times[i_credit]
        qs = curve_qs[i_credit]

        for i_trial in range(num_trials):

            corr_times[i_credit, i_trial] = (
                uniform_to_default_time(
                    u1[i_credit, i_trial], times, qs
                )
            )

            corr_times[i_credit, i_trial + num_trials] = (
                uniform_to_default_time(
                    u2[i_credit, i_trial], times, qs
                )
            )

    return corr_times


@njit(fastmath=True, cache=True)
def uniform_to_default_time(u, t, v):
    """Fast mapping of a uniform random variable to a default time given a
    survival probability curve."""

    if u <= 0.0:
        return 99999.0

    if u >= 1.0:
        return 0.0

    num_points = len(v)
    index = 0

    for i in range(1, num_points):
        if u <= v[i - 1] and u > v[i]:
            index = i
            break

    if index == 0:
        # Extrapolate beyond the final curve node using its last forward hazard.
        left = num_points - 2
        right = num_points - 1
        lam = np.log(v[left] / v[right]) / (t[right] - t[left])
        if lam <= 0.0:
            return 99999.0
        tau = t[right] + np.log(v[right] / u) / lam
    else:
        t1 = t[index - 1]
        q1 = v[index - 1]
        t2 = t[index]
        q2 = v[index]
        lam = np.log(q1 / q2) / (t2 - t1)
        tau = t1 + np.log(q1 / u) / lam

    return tau
