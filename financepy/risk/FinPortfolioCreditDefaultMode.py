##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp
import numpy as np
from ..finutils.FinHelperFunctions import normaliseWeights
from ..models.FinModelGaussianCopula1F import lossDbnHeterogeneousAdjBinomial


class FinPortfolioCreditDefaultMode(object):

    def __init__(self,
                 weights):

        self._numCredits = len(weights)
        self._weights = weights

###############################################################################

    def lossDistribution(self,
                         tmat,
                         hazardRates,
                         recoveryRates,
                         betaValues,
                         numPoints):

        if tmat < 0.0:
            raise ValueError("Maturity must be positive")

        if len(hazardRates) != self._numCredits:
            raise ValueError("Hazard rates dop not have " /
                             + str(self._numCredits) + "elements.")

        if len(recoveryRates) != self._numCredits:
            raise ValueError("Recovery rates dop not have " /
                             + str(self._numCredits) + "elements.")

        if len(betaValues) != self._numCredits:
            raise ValueError("Beta values do not have " /
                             + str(self._numCredits) + "elements.")

        if numPoints < 1:
            raise ValueError("Num points must be > 0")

        self._hazardRates = hazardRates
        self._recoveryRates = recoveryRates
        self._betaValues = betaValues

        weights = normaliseWeights(self._weights)
        defaultProbs = np.zeros(self._numCredits)

        for j in range(0, self._numCredits):
            defaultProbs[j] = 1.0 - exp(-hazardRates[j] * tmat)

        totalLoss = 0.0
        for i in range(0, self._numCredits):
            totalLoss = totalLoss + weights[i] * (1.0 - recoveryRates[i])
        avgLoss = totalLoss / self._numCredits

        lossRatio = np.zeros(self._numCredits)
        for i in range(0, self._numCredits):
            lossRatio[i] = weights[i] * (1.0 - recoveryRates[i]) / avgLoss

        self._support = np.linspace(0.0, totalLoss, self._numCredits+1)
        self._lossDbn = lossDbnHeterogeneousAdjBinomial(self._numCredits,
                                                        defaultProbs,
                                                        lossRatio,
                                                        betaValues,
                                                        numPoints)

        return (self._support, self._lossDbn)

###############################################################################
