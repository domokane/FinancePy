##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

from ..utils.math import pair_gcd

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_heterogeneous_adj_binomial(num_credits,
                                              condProbs,
                                              loss_ratio):

    # Algorithm due to D. O'Kane.

    numLosses = num_credits + 1
    indepDbn = np.zeros(numLosses)

    p = 0.0
    for iCredit in range(0, num_credits):
        p += loss_ratio[iCredit] * condProbs[iCredit]
    p = p / num_credits

    ###########################################################################

    if p < 0.5:
        ratio = p / (1.0 - p)
        indepDbn[0] = (1.0 - p)**num_credits
        for i in range(1, numLosses):
            indepDbn[i] = indepDbn[i - 1] * ratio * (num_credits - i + 1.0) / i
    else:
        ratio = (1.0 - p) / p
        indepDbn[num_credits] = p ** num_credits
        for i in range(num_credits - 1, -1, -1):
            indepDbn[i] = indepDbn[i + 1] * \
                ratio * (i + 1.0) / (num_credits - i)

    ###########################################################################

    vapprox = 0.0
    vexact = 0.0

    for iCredit in range(0, num_credits):
        loss_ratio2 = loss_ratio[iCredit] ** 2
        vapprox += loss_ratio2 * p * (1.0 - p)
        vexact += loss_ratio2 * condProbs[iCredit] * (1.0 - condProbs[iCredit])

    ###########################################################################

    meanLoss = p * num_credits
    meanAbove = round(meanLoss + 1)
    meanBelow = round(meanLoss)

    if meanAbove > num_credits:
        meanAbove = num_credits

    diffAbove = meanAbove - meanLoss
    diffBelow = meanBelow - meanLoss

    # DOK - TO DO - SIMPLIFY THIS CODE AS PER JOD PAPER
    term = diffAbove * diffAbove + \
        (diffBelow * diffBelow - diffAbove * diffAbove) * diffAbove
    numer = vexact - term
    denom = vapprox - term

    if abs(denom) < 1e-30:
        denom = 1e-30

    alpha = numer / denom
    epsilonBelow = (1.0 - alpha) * diffAbove
    epsilonAbove = (1.0 - alpha) - epsilonBelow

    for iLossUnit in range(0, numLosses):
        indepDbn[iLossUnit] *= alpha

    indepDbn[int(meanBelow)] += epsilonBelow
    indepDbn[int(meanAbove)] += epsilonAbove

    return indepDbn

###############################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def portfolio_gcd(actualLosses):

    num_credits = len(actualLosses)
    scaling = 1000000

    temp = (int)(actualLosses[0] * scaling)

    for iCredit in range(1, num_credits):
        num2 = int(actualLosses[iCredit] * scaling)
        temp = pair_gcd(temp, num2)

    portfolioGCD = float(temp / scaling)
    return portfolioGCD

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_recursion_gcd(num_credits,
                                 condDefaultProbs,
                                 lossUnits):

    numLossUnits = 1
    for i in range(0, len(lossUnits)):
        numLossUnits += int(lossUnits[i])

    prevDbn = np.zeros(numLossUnits)
    prevDbn[0] = 1.0

    small = 1e-10
    nextDbn = np.zeros(numLossUnits)

    for iCredit in range(0, num_credits):

        p = condDefaultProbs[iCredit]
        loss = (int)(lossUnits[iCredit] + small)

        for iLossUnit in range(0, loss):
            nextDbn[iLossUnit] = prevDbn[iLossUnit] * (1.0 - p)

        for iLossUnit in range(loss, numLossUnits):
            nextDbn[iLossUnit] = prevDbn[iLossUnit - loss] * \
                p + prevDbn[iLossUnit] * (1.0 - p)

        for iLossUnit in range(0, numLossUnits):
            prevDbn[iLossUnit] = nextDbn[iLossUnit]

    return nextDbn

##########################################################################
