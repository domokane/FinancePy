##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt
import numpy as np
from scipy.stats import t as student

from ..finutils.FinHelperFunctions import uniformToDefaultTime

###############################################################################


class FinModelStudentTCopula():

    def defaultTimes(self,
                     issuerCurves,
                     correlationMatrix,
                     degreesOfFreedom,
                     numTrials,
                     seed):

        np.random.seed(seed)
        numCredits = len(issuerCurves)
        x = np.random.normal(0.0, 1.0, size=(numCredits, numTrials))
        c = np.linalg.cholesky(correlationMatrix)
        y = np.dot(c, x)

        corrTimes = np.empty(shape=(numCredits, 2 * numTrials))

        for iTrial in range(0, numTrials):
            chi2 = np.random.chisquare(degreesOfFreedom)
            c = sqrt(chi2 / degreesOfFreedom)
            for iCredit in range(0, numCredits):
                issuerCurve = issuerCurves[iCredit]
                g = y[iCredit, iTrial] / c
                u1 = student.cdf(g, degreesOfFreedom)
                u2 = 1.0 - u1
                times = issuerCurve._times
                values = issuerCurve._values
                t1 = uniformToDefaultTime(u1, times, values)
                t2 = uniformToDefaultTime(u2, times, values)
                corrTimes[iCredit, iTrial] = t1
                corrTimes[iCredit, iTrial + numTrials] = t2

        return corrTimes

###############################################################################
