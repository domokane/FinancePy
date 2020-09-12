##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Use Numba to speed this up more.

import numpy as np

from ..finutils.FinMath import N
from ..finutils.FinHelperFunctions import uniformToDefaultTime

###############################################################################
# TODO:
###############################################################################


def defaultTimesGC(issuerCurves,
                   correlationMatrix,
                   numTrials,
                   seed):
    ''' Generate a matrix of default times by credit and trial using a
    Gaussian copula model using a full rank correlation matrix. '''

    np.random.seed(seed)
    numCredits = len(issuerCurves)
    x = np.random.normal(0.0, 1.0, size=(numCredits, numTrials))
    c = np.linalg.cholesky(correlationMatrix)
    y = np.dot(c, x)

    corrTimes = np.empty(shape=(numCredits, 2 * numTrials))

    for iCredit in range(0, numCredits):
        issuerCurve = issuerCurves[iCredit]
        for iTrial in range(0, numTrials):
            g = y[iCredit, iTrial]
            u1 = 1.0 - N(g)
            u2 = 1.0 - u1
            times = issuerCurve._times
            values = issuerCurve._values
            t1 = uniformToDefaultTime(u1, times, values)
            t2 = uniformToDefaultTime(u2, times, values)
            corrTimes[iCredit, iTrial] = t1
            corrTimes[iCredit, numTrials + iTrial] = t2

    return corrTimes

##########################################################################
