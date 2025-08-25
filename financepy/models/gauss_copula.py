# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# TODO: Use Numba to speed this up more.

import numpy as np

from ..utils.math import normcdf
from ..utils.helpers import uniform_to_default_time

# TODO:

########################################################################################


def default_times_gc(issuer_curves, corr_matrix, num_trials, seed):
    """Generate a matrix of default times by credit and trial using a
    Gaussian copula model using a full rank correlation matrix."""

    np.random.seed(seed)
    num_credits = len(issuer_curves)
    x = np.random.normal(0.0, 1.0, size=(num_credits, num_trials))
    c = np.linalg.cholesky(corr_matrix)
    y = np.dot(c, x)

    corr_times = np.empty(shape=(num_credits, 2 * num_trials))

    for i_credit in range(0, num_credits):
        issuer_curve = issuer_curves[i_credit]
        for i_trial in range(0, num_trials):
            g = y[i_credit, i_trial]
            u1 = 1.0 - normcdf(g)
            u2 = 1.0 - u1
            times = issuer_curve.times
            qs = issuer_curve.qs
            t1 = uniform_to_default_time(u1, times, qs)
            t2 = uniform_to_default_time(u2, times, qs)
            corr_times[i_credit, i_trial] = t1
            corr_times[i_credit, num_trials + i_trial] = t2

    return corr_times
