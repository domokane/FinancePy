# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 21:57:13 2019

@author: Dominic
"""

from numba import njit, float64

###############################################################################
# TODO: Incorporate interpolation at boundary
###############################################################################


@njit(float64(float64[:], float64[:]), fastmath=True, cache=True)
def expectedLoss(lossSizeVector,
                 lossProbabilityVector):

    numLossUnits = lossProbabilityVector.shape[0]

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if sumProbs != 1.0:
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    expectedLoss = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        loss = lossSizeVector[i]
        expectedLoss = expectedLoss + prob * loss

    return expectedLoss

###############################################################################


@njit(fastmath=True, cache=True)
def valueAtRisk(lossSizeVector,
                lossProbabilityVector,
                confidenceLevel):

    if len(lossSizeVector) < 2:
        raise ValueError("Loss vector requires 2 or more entries.")

    if confidenceLevel < 0.0 or confidenceLevel > 1.0:
        raise ValueError("Confidence level must be [0,1]")

    numLossUnits = len(lossProbabilityVector)

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if sumProbs != 1.0:
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    cumProb = 0.0
    iExceed = 0
    for i in range(0, numLossUnits):
        cumProb += lossProbabilityVector[i]
        if cumProb > confidenceLevel:
            iExceed = i
            break

    # Interpolation of the VaR based on the cumulative probability
    cumProbBelow = cumProb - lossProbabilityVector[iExceed]
    var = (1.0 - confidenceLevel) * lossSizeVector[iExceed-1]
    var = var + (confidenceLevel - cumProbBelow) * lossSizeVector[iExceed]
    var = var / (1.0 - cumProbBelow)
    return var

###############################################################################


@njit(fastmath=True, cache=True)
def expectedShortfall(lossSizeVector,
                      lossProbabilityVector,
                      confidenceLevel):

    if confidenceLevel < 0.0 or confidenceLevel > 1.0:
        raise ValueError("Confidence level must be [0,1]")

    if len(lossSizeVector) < 2:
        raise ValueError("Loss vector requires 2 or more entries.")

    numLossUnits = len(lossProbabilityVector)

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if sumProbs != 1.0:
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    cumProb = 0.0
    iExceed = 0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        cumProb += prob
        if cumProb > confidenceLevel:
            iExceed = i
            break

    # TO DO - GET THE INTERPOLATION WORKING HERE !
    expectedShortfall = 0.0
    tailProbability = 0.0
    for i in range(iExceed, numLossUnits):
        prob = lossProbabilityVector[i]
        loss = lossSizeVector[i]
        tailProbability += prob
        expectedShortfall += prob * loss

    if tailProbability > 0.0:
        expectedShortfall /= tailProbability
    else:
        raise ValueError("Tail probability is zero. How is this possible ?")

    return expectedShortfall

###############################################################################
