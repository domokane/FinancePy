##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

##########################################################################

from ..finutils.FinMath import norminvcdf, N, INVROOT2PI
from ..finutils.FinError import FinError
from .FinModelLossDbnBuilder import indepLossDbnRecursionGCD
from .FinModelLossDbnBuilder import indepLossDbnHeterogeneousAdjBinomial
from .FinModelLossDbnBuilder import portfolioGCD

###############################################################################

minZ = -6.0

###############################################################################
# This implements the one-factor latent variable formulation of the Gaussian
# Copula model as well as some approximations
###############################################################################

@njit(float64[:](int64, float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def lossDbnRecursionGCD(numCredits,
                        defaultProbs,
                        lossUnits,
                        betaVector,
                        numIntegrationSteps):
    ''' Full construction of the loss distribution of a portfolio of credits
    where losses have been calculate as number of units based on the GCD. '''

    if len(defaultProbs) != numCredits:
        raise FinError("Default probability length must equal num credits.")

    if len(lossUnits) != numCredits:
        raise FinError("Loss units length must equal num credits.")

    if len(betaVector) != numCredits:
        raise FinError("Beta vector length must equal num credits.")

    numLossUnits = 1
    for i in range(0, len(lossUnits)):
        numLossUnits += int(lossUnits[i])

    uncondLossDbn = np.zeros(numLossUnits)

    z = minZ
    dz = 2.0 * abs(z) / numIntegrationSteps

    condDefaultProbs = np.zeros(numCredits)

    thresholds = np.zeros(numCredits)
    for iCredit in range(0, int(numCredits)):
        thresholds[iCredit] = norminvcdf(defaultProbs[iCredit])

    for _ in range(0, numIntegrationSteps):
        for iCredit in range(0, numCredits):
            beta = betaVector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indepLossDbnRecursionGCD(numCredits,
                                            condDefaultProbs,
                                            lossUnits)

        gaussWt = np.exp(-(z*z)/2.0)

        for iLossUnit in range(0, numLossUnits):
            uncondLossDbn[iLossUnit] += indepDbn[iLossUnit] * gaussWt

        z += dz

    for iLossUnit in range(0, int(numLossUnits)):
        uncondLossDbn[iLossUnit] *= INVROOT2PI * dz

    return uncondLossDbn

###############################################################################

@njit(float64[:](float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def homogeneousBasketLossDbn(survivalProbabilities,
                             recoveryRates,
                             betaVector,
                             numIntegrationSteps):
    ''' Calculate the loss distribution of a CDS default basket where the
    portfolio is equally weighted and the losses in the portfolio are homo-
    geneous i.e. the credits have the same recovery rates. '''

    numCredits = len(survivalProbabilities)

    if numCredits == 0:
        raise FinError("Number of credits equals zero")

    for iCredit in range(1, numCredits):
        if recoveryRates[iCredit] != recoveryRates[0]:
            raise FinError("Losses are not homogeneous")

    m = 0.0
    for i in range(0, len(betaVector)):
        m += betaVector[i]
    m /= len(betaVector)

    # High beta requires more integration steps
    if m > 0.7:
        numIntegrationSteps *= 2

    if m > 0.9:
        numIntegrationSteps *= 5

    lossUnits = np.ones(numCredits)

    defaultProbs = np.zeros(numCredits)
    for iCredit in range(0, numCredits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    lossDbn = lossDbnRecursionGCD(numCredits,
                                  defaultProbs,
                                  lossUnits,
                                  betaVector,
                                  numIntegrationSteps)

    return lossDbn

###############################################################################


@njit(float64(float64, float64, int64, float64[:], float64[:], float64[:],
              int64), fastmath=True)
def trSurvProbRecursion(k1,
                        k2,
                        numCredits,
                        survivalProbabilities,
                        recoveryRates,
                        betaVector,
                        numIntegrationSteps):
    ''' Get the tranche survival probability of a portfolio of credits in the
    one-factor GC model using a full recursion calculation of the loss
    distribution and survival probabilities to some time horizon. '''

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    commonRecoveryFlag = 1

    lossAmounts = np.zeros(numCredits)
    for iCredit in range(0, numCredits):
        lossAmounts[iCredit] = (1.0 - recoveryRates[iCredit]) / numCredits
        if lossAmounts[iCredit] != lossAmounts[0]:
            commonRecoveryFlag = 0

    gcd = 0.0

    m = 0.0
    for i in range(0, len(betaVector)):
        m += betaVector[i]
    m /= len(betaVector)

    if m > 0.8:
        numIntegrationSteps *= 2

    if commonRecoveryFlag == 1:
        gcd = lossAmounts[0]
    else:
        gcd = portfolioGCD(lossAmounts)

    lossUnits = np.zeros(numCredits)
    numLossUnits = 1  # this is the zero loss

    for iCredit in range(0, numCredits):
        lossUnits[iCredit] = lossAmounts[iCredit] / gcd
        numLossUnits = numLossUnits + lossUnits[iCredit]

    defaultProbs = np.zeros(numCredits)

    for iCredit in range(0, numCredits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    lossDbn = lossDbnRecursionGCD(numCredits,
                                  defaultProbs,
                                  lossUnits,
                                  betaVector,
                                  numIntegrationSteps)

    trancheEL = 0.0
    for iLossUnit in range(0, int(numLossUnits)):
        loss = iLossUnit * gcd
        trancheLoss = min(loss, k2) - min(loss, k1)
        trancheEL = trancheEL + trancheLoss * lossDbn[iLossUnit]

    q = 1.0 - trancheEL / (k2 - k1)
    return q

###############################################################################


@njit(float64(float64, float64, float64, float64), fastmath=True, cache=True)
def gaussApproxTrancheLoss(k1, k2, mu, sigma):

    if abs(sigma) < 1e-6:
        gaussApproxTrancheLoss = 0.0
        if mu > k1:
            gaussApproxTrancheLoss += (mu - k1)

        if mu > k2:
            gaussApproxTrancheLoss += (mu - k2)
    else:
        d1 = (mu - k1) / sigma
        d2 = (mu - k2) / sigma

        gaussApproxTrancheLoss = (mu - k1) * N(d1) - (mu - k2) * N(d2)
        + sigma * np.exp(-0.5 * d1 * d1) * INVROOT2PI
        - sigma * np.exp(-0.5 * d2 * d2) * INVROOT2PI

    return gaussApproxTrancheLoss

###############################################################################


@njit(float64(float64, float64, int64, float64[:], float64[:], float64[:],
              int64), fastmath=True, cache=True)
def trSurvProbGaussian(k1,
                       k2,
                       numCredits,
                       survivalProbabilities,
                       recoveryRates,
                       betaVector,
                       numIntegrationSteps):
    ''' Get the approximated tranche survival probability of a portfolio
    of credits in the one-factor GC model using a Gaussian fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. Note that the losses in this fit are allowed to be negative. '''

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    defaultProbs = [0.0] * numCredits
    for iCredit in range(0, numCredits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    dz = 2.0 * abs(minZ) / numIntegrationSteps
    z = minZ

    thresholds = np.zeros(numCredits)
    losses = np.zeros(numCredits)

    for iCredit in range(0, numCredits):
        pd = 1.0 - survivalProbabilities[iCredit]
        thresholds[iCredit] = norminvcdf(pd)
        losses[iCredit] = (1.0 - recoveryRates[iCredit]) / numCredits

    v = 0.0
    for _ in range(0, numIntegrationSteps):

        mu = 0.0
        var = 0.0

        # calculate the mean and variance of the conditional loss distribution
        for iCredit in range(0, numCredits):
            beta = betaVector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condprob = N(argz)
            mu += condprob * losses[iCredit]
            var += (losses[iCredit]**2) * condprob * (1.0 - condprob)

        sigma = np.sqrt(var)
        el = gaussApproxTrancheLoss(k1, k2, mu, sigma)
        gaussWt = np.exp(-(z**2) / 2.0)

        v += el * gaussWt
        z += dz

    v *= INVROOT2PI * dz
    q = 1.0 - v / (k2 - k1)
    return q

###############################################################################

@njit(float64[:](int64, float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def lossDbnHeterogeneousAdjBinomial(numCredits,
                                    defaultProbs,
                                    lossRatio,
                                    betaVector,
                                    numIntegrationSteps):
    ''' Get the portfolio loss distribution using the adjusted binomial
    approximation to the conditional loss distribution. '''

    numLossUnits = numCredits + 1
    condDefaultProbs = np.zeros(numCredits)
    uncondLossDbn = np.zeros(numLossUnits)
    indepDbn = np.zeros(numLossUnits)

    # Determine default threshold for each credit
    thresholds = np.zeros(numCredits)
    for iCredit in range(0, numCredits):
        thresholds[iCredit] = norminvcdf(defaultProbs[iCredit])

    dz = 2.0 * abs(minZ) / numIntegrationSteps
    z = minZ

    for iStep in range(0, numIntegrationSteps):
        for iCredit in range(0, numCredits):
            beta = betaVector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indepLossDbnHeterogeneousAdjBinomial(numCredits,
                                                        condDefaultProbs,
                                                        lossRatio)

        gaussWt = np.exp(-(z**2) / 2.0)

        for iLossUnit in range(0, numLossUnits):
            uncondLossDbn[iLossUnit] += indepDbn[iLossUnit] * gaussWt

        z = z + dz

    for iLossUnit in range(0, numLossUnits):
        uncondLossDbn[iLossUnit] *= INVROOT2PI * dz

    return uncondLossDbn

###############################################################################


@njit(float64(float64, float64, int64, float64[:], float64[:], float64[:],
              int64), fastmath=True, cache=True)
def trSurvProbAdjBinomial(k1,
                          k2,
                          numCredits,
                          survivalProbabilities,
                          recoveryRates,
                          betaVector,
                          numIntegrationSteps):
    ''' Get the approximated tranche survival probability of a portfolio of
    credits in the one-factor GC model using the adjusted binomial fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. This approach is both fast and highly accurate. '''

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    defaultProbs = np.zeros(numCredits)
    for iCredit in range(0, numCredits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    totalLoss = 0.0
    for iCredit in range(0, numCredits):
        totalLoss += (1.0 - recoveryRates[iCredit])
    totalLoss /= numCredits

    avgLoss = totalLoss / numCredits

    lossRatio = np.zeros(numCredits)
    for iCredit in range(0, numCredits):
        lossRatio[iCredit] = (
            1.0 - recoveryRates[iCredit]) / numCredits / avgLoss

    lossDbn = lossDbnHeterogeneousAdjBinomial(numCredits,
                                              defaultProbs,
                                              lossRatio,
                                              betaVector,
                                              numIntegrationSteps)
    trancheEL = 0.0
    numLossUnits = numCredits + 1
    for iLossUnit in range(0, numLossUnits):
        loss = iLossUnit * avgLoss
        trancheLoss = min(loss, k2) - min(loss, k1)
        trancheEL += trancheLoss * lossDbn[iLossUnit]

    q = 1.0 - trancheEL / (k2 - k1)
    return q

###############################################################################
