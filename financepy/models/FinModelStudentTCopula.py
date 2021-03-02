##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt
import numpy as np
from scipy.stats import t as student

from ..utils.FinHelperFunctions import uniformToDefaultTime

###############################################################################


class FinModelStudentTCopula():

    def defaultTimes(self,
                     issuer_curves,
                     correlationMatrix,
                     degreesOfFreedom,
                     numTrials,
                     seed):

        np.random.seed(seed)
        numCredits = len(issuer_curves)
        x = np.random.normal(0.0, 1.0, size=(numCredits, numTrials))
        c = np.linalg.cholesky(correlationMatrix)
        y = np.dot(c, x)

        corrTimes = np.empty(shape=(numCredits, 2 * numTrials))

        for iTrial in range(0, numTrials):
            chi2 = np.random.chisquare(degreesOfFreedom)
            c = sqrt(chi2 / degreesOfFreedom)
            for iCredit in range(0, numCredits):
                issuer_curve = issuer_curves[iCredit]
                g = y[iCredit, iTrial] / c
                u1 = student.cdf(g, degreesOfFreedom)
                u2 = 1.0 - u1
                times = issuer_curve._times
                values = issuer_curve._values
                t1 = uniformToDefaultTime(u1, times, values)
                t2 = uniformToDefaultTime(u2, times, values)
                corrTimes[iCredit, iTrial] = t1
                corrTimes[iCredit, iTrial + numTrials] = t2

        return corrTimes

###############################################################################
