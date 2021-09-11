##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Use Numba to speed this up more.

import numpy as np

from ..utils.math import N
from ..utils.helpers import uniform_to_default_time

###############################################################################
# TODO:
###############################################################################


def default_times_gc(issuer_curves,
                     correlationMatrix,
                     num_trials,
                     seed):
    """ Generate a matrix of default times by credit and trial using a
    Gaussian copula model using a full rank correlation matrix. """

    np.random.seed(seed)
    num_credits = len(issuer_curves)
    x = np.random.normal(0.0, 1.0, size=(num_credits, num_trials))
    c = np.linalg.cholesky(correlationMatrix)
    y = np.dot(c, x)

    corrTimes = np.empty(shape=(num_credits, 2 * num_trials))

    for iCredit in range(0, num_credits):
        issuer_curve = issuer_curves[iCredit]
        for iTrial in range(0, num_trials):
            g = y[iCredit, iTrial]
            u1 = 1.0 - N(g)
            u2 = 1.0 - u1
            times = issuer_curve._times
            values = issuer_curve._values
            t1 = uniform_to_default_time(u1, times, values)
            t2 = uniform_to_default_time(u2, times, values)
            corrTimes[iCredit, iTrial] = t1
            corrTimes[iCredit, num_trials + iTrial] = t2

    return corrTimes

##########################################################################
