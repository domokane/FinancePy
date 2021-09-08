##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt
import numpy as np
from scipy.stats import t as student

from ..utils.helpers import uniform_to_default_time

###############################################################################


class StudentTCopula():

    def default_times(self,
                      issuer_curves,
                      correlationMatrix,
                      degreesOfFreedom,
                      num_trials,
                      seed):

        np.random.seed(seed)
        num_credits = len(issuer_curves)
        x = np.random.normal(0.0, 1.0, size=(num_credits, num_trials))
        c = np.linalg.cholesky(correlationMatrix)
        y = np.dot(c, x)

        corrTimes = np.empty(shape=(num_credits, 2 * num_trials))

        for iTrial in range(0, num_trials):
            chi2 = np.random.chisquare(degreesOfFreedom)
            c = sqrt(chi2 / degreesOfFreedom)
            for iCredit in range(0, num_credits):
                issuer_curve = issuer_curves[iCredit]
                g = y[iCredit, iTrial] / c
                u1 = student.cdf(g, degreesOfFreedom)
                u2 = 1.0 - u1
                times = issuer_curve._times
                values = issuer_curve._values
                t1 = uniform_to_default_time(u1, times, values)
                t2 = uniform_to_default_time(u2, times, values)
                corrTimes[iCredit, iTrial] = t1
                corrTimes[iCredit, iTrial + num_trials] = t2

        return corrTimes

###############################################################################
