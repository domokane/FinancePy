# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# TODO: Use Numba to speed this up more.

from scipy.stats import qmc
from scipy.special import ndtri
from scipy.special import ndtr

from numba.typed import List
import math
import numpy as np
from .default_times import default_times_from_uniforms_numba

from typing import Sequence, Any


def gauss_default_times(
        issuer_curves: Sequence[Any],
        corr_matrix: np.ndarray,
        num_trials: int,
        seed: int,
        use_moment_matching: bool = True,
        use_sobol: bool = False):

    num_credits = len(issuer_curves)

    if use_sobol:
        sampler = qmc.Sobol(
            d=num_credits,
            scramble=True,
            seed=seed
        )

        m = (num_trials - 1).bit_length()
        u = sampler.random_base2(m)

        num_trials = len(u)

        x = ndtri(np.clip(u, 1e-15, 1.0 - 1e-15)).T

    else:

        rng = np.random.RandomState(seed)

        x = rng.normal(
            0.0, 1.0, size=(num_credits, num_trials)
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

    u1 = ndtr(-y)
    u2 = ndtr(y)

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
