##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from numba import njit, float64
from ..finutils.FinMath import testMonotonicity

small = 1e-8
###############################################################################
# TODO: Incorporate interpolation at boundary
###############################################################################


@njit(float64(float64[:], float64[:]), fastmath=True, cache=True)
def expectedLoss(lossSizeVector,
                 lossProbabilityVector):

    numLossUnits = len(lossProbabilityVector)

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    if testMonotonicity(lossSizeVector) is False:
        raise ValueError("Losses are not monotonic.")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if abs(sumProbs - 1.0) > small:
        print(sumProbs)
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    expectedLoss = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        loss = lossSizeVector[i]
        expectedLoss = expectedLoss + prob * loss

    return expectedLoss

###############################################################################


@njit(float64(float64[:], float64[:], float64),fastmath=True, cache=True)
def valueAtRisk(lossSizeVector,
                lossProbabilityVector,
                confidenceLevel):

    if len(lossSizeVector) < 2:
        raise ValueError("Loss vector requires 2 or more entries.")

    if confidenceLevel < 0.0 or confidenceLevel > 1.0:
        raise ValueError("Confidence level must be [0,1]")

    if testMonotonicity(lossSizeVector) is False:
        raise ValueError("Losses are not monotonic.")

    numLossUnits = len(lossProbabilityVector)

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if abs(sumProbs - 1.0) > small:
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    cumProb = 0.0
    iExceed = -1
    for i in range(0, numLossUnits):
        cumProb += lossProbabilityVector[i]
        if cumProb > confidenceLevel:
            iExceed = i
            break

    # If the VaR is in the first spike of the loss distribution
    if iExceed == 0:
        return lossSizeVector[0]

    # if the VaR is in the last spike of the loss distribution then find loss
    if iExceed == -1:
        iMaxLoss = 0
        for i in range(0, numLossUnits):
            if lossProbabilityVector[i] > 0.0:
                iMaxLoss = i
        return lossSizeVector[iMaxLoss]

    var = lossSizeVector[iExceed]
    return var

###############################################################################


@njit(float64(float64[:], float64[:], float64), fastmath=True, cache=True)
def expectedShortfall(lossSizeVector,
                      lossProbabilityVector,
                      confidenceLevel):

    if confidenceLevel < 0.0 or confidenceLevel > 1.0:
        raise ValueError("Confidence level must be [0,1]")

    if len(lossSizeVector) < 2:
        raise ValueError("Loss vector requires 2 or more entries.")

    if testMonotonicity(lossSizeVector) is False:
        raise ValueError("Losses are not monotonic.")

    numLossUnits = len(lossProbabilityVector)

    if len(lossSizeVector) != numLossUnits:
        raise ValueError("Loss amount and probability vectors not same length")

    sumProbs = 0.0
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        sumProbs += prob

    if abs(sumProbs - 1.0) > small:
        raise ValueError("Loss distribution probabilities do not sum to 1.0")

    cumProb = 0.0
    iExceed = -1
    for i in range(0, numLossUnits):
        prob = lossProbabilityVector[i]
        cumProb += prob
        if cumProb > confidenceLevel:
            iExceed = i
            break

    # We split the first spike so that we only count bit inside tail
    expectedLoss = (cumProb - confidenceLevel) * lossSizeVector[iExceed]

    for i in range(iExceed+1, numLossUnits):
        prob = lossProbabilityVector[i]
        loss = lossSizeVector[i]
        expectedLoss += prob * loss

    tailProbability = 1.0 - confidenceLevel

    if tailProbability > 0.0:
        expectedShortfall = expectedLoss / tailProbability
    else:
        expectedShortfall = lossSizeVector[-1]

    return expectedShortfall

###############################################################################
