# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from math import sqrt
import numpy as np
from scipy.stats import t as student

from ..utils.helpers import uniform_to_default_time

########################################################################################



from typing import Sequence, Any

class StudentTCopula:

    ####################################################################################

    def default_times(
        self,
        issuer_curves: Sequence[Any],
        corr_matrix: np.ndarray,
        degrees_of_freedom: float,
        num_trials: int,
        seed: int
    ) -> np.ndarray:

        np.random.seed(seed)
        num_credits = len(issuer_curves)
        x = np.random.normal(0.0, 1.0, size=(num_credits, num_trials))
        c = np.linalg.cholesky(corr_matrix)
        y = np.dot(c, x)

        corr_times = np.empty(shape=(num_credits, 2 * num_trials))

        for i_trial in range(0, num_trials):
            chi2 = np.random.chisquare(degrees_of_freedom)
            c = sqrt(chi2 / degrees_of_freedom)
            for i_credit in range(0, num_credits):
                issuer_curve = issuer_curves[i_credit]
                g = y[i_credit, i_trial] / c
                u1 = student.cdf(g, degrees_of_freedom)
                u2 = 1.0 - u1
                times = issuer_curve.times
                qs = issuer_curve.qs
                t1 = uniform_to_default_time(u1, times, qs)
                t2 = uniform_to_default_time(u2, times, qs)
                corr_times[i_credit, i_trial] = t1
                corr_times[i_credit, i_trial + num_trials] = t2

        return corr_times
