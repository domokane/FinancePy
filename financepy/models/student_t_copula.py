# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from scipy.stats import qmc
from scipy.special import ndtri
from scipy.stats import chi2

import numpy as np
from scipy.stats import t as student
from numba.typed import List
from .default_times import default_times_from_uniforms_numba


from typing import Sequence, Any

########################################################################################


def student_t_default_times(
        issuer_curves: Sequence[Any],
        corr_matrix: np.ndarray,
        degrees_of_freedom: float,
        num_trials: int,
        seed: int,
        use_moment_matching: bool = True,
        use_sobol: bool = False):

    num_credits = len(issuer_curves)

    rng = np.random.RandomState(seed)

    if use_sobol:
        sampler = qmc.Sobol(
            d=num_credits+1,
            scramble=True,
            seed=seed
        )

        m = (num_trials - 1).bit_length()
        u = sampler.random_base2(m)

        num_trials = len(u)

        u = np.clip(u, 1e-15, 1.0 - 1e-15)

        x = ndtri(u[:, :num_credits]).T

        chi2_draws = chi2.ppf(
            u[:, num_credits],
            degrees_of_freedom
        )

    else:

        x = rng.normal(
            0.0, 1.0, size=(num_credits, num_trials)
        )

        chi2_draws = rng.chisquare(
            degrees_of_freedom,
            size=num_trials
        )

    # Moment matching
    if use_moment_matching:
        x -= np.mean(x, axis=1, keepdims=True)
        x /= np.std(x, axis=1, keepdims=True)

    if np.allclose(corr_matrix, 1.0, rtol=0.0, atol=1e-14):
        # Perfect positive correlation: one common Gaussian factor
        y = np.broadcast_to(x[0:1, :], x.shape).copy()
    else:
        c = np.linalg.cholesky(corr_matrix)
        y = c @ x

    # Apply the common chi-squared scaling
    c = np.sqrt(chi2_draws / degrees_of_freedom)

    g = y / c[np.newaxis, :]

    # Vectorised Student's t CDF
    u1 = student.cdf(g, degrees_of_freedom)
    u2 = 1.0 - u1

    # Invert survival curves using Numba
    curve_times = List()
    curve_qs = List()

    for curve in issuer_curves:
        curve_times.append(
            np.asarray(curve._times, dtype=np.float64)
        )
        curve_qs.append(
            np.asarray(curve._qs, dtype=np.float64)
        )

    return default_times_from_uniforms_numba(
        u1, u2, curve_times, curve_qs
    )
