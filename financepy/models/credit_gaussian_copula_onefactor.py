##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

##########################################################################

from ..utils.math import norminvcdf, N, INVROOT2PI
from ..utils.error import FinError
from .credit_loss_dbn_builder import indepLossDbnRecursionGCD
from .credit_loss_dbn_builder import indepLossDbnHeterogeneousAdjBinomial
from .credit_loss_dbn_builder import portfolioGCD

###############################################################################

minZ = -6.0

###############################################################################
# This implements the one-factor latent variable formulation of the Gaussian
# Copula model as well as some approximations
###############################################################################

@njit(float64[:](int64, float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def lossDbnRecursionGCD(num_credits,
                        defaultProbs,
                        lossUnits,
                        betaVector,
                        numIntegrationSteps):
    """ Full construction of the loss distribution of a portfolio of credits
    where losses have been calculate as number of units based on the GCD. """

    if len(defaultProbs) != num_credits:
        raise FinError("Default probability length must equal num credits.")

    if len(lossUnits) != num_credits:
        raise FinError("Loss units length must equal num credits.")

    if len(betaVector) != num_credits:
        raise FinError("Beta vector length must equal num credits.")

    numLossUnits = 1
    for i in range(0, len(lossUnits)):
        numLossUnits += int(lossUnits[i])

    uncondLossDbn = np.zeros(numLossUnits)

    z = minZ
    dz = 2.0 * abs(z) / numIntegrationSteps

    condDefaultProbs = np.zeros(num_credits)

    thresholds = np.zeros(num_credits)
    for iCredit in range(0, int(num_credits)):
        thresholds[iCredit] = norminvcdf(defaultProbs[iCredit])

    for _ in range(0, numIntegrationSteps):
        for iCredit in range(0, num_credits):
            beta = betaVector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indepLossDbnRecursionGCD(num_credits,
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
                             recovery_rates,
                             betaVector,
                             numIntegrationSteps):
    """ Calculate the loss distribution of a CDS default basket where the
    portfolio is equally weighted and the losses in the portfolio are homo-
    geneous i.e. the credits have the same recovery rates. """

    num_credits = len(survivalProbabilities)

    if num_credits == 0:
        raise FinError("Number of credits equals zero")

    for iCredit in range(1, num_credits):
        if recovery_rates[iCredit] != recovery_rates[0]:
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

    lossUnits = np.ones(num_credits)

    defaultProbs = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    lossDbn = lossDbnRecursionGCD(num_credits,
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
                        num_credits,
                        survivalProbabilities,
                        recovery_rates,
                        betaVector,
                        numIntegrationSteps):
    """ Get the tranche survival probability of a portfolio of credits in the
    one-factor GC model using a full recursion calculation of the loss
    distribution and survival probabilities to some time horizon. """

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    commonRecoveryFlag = 1

    lossAmounts = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        lossAmounts[iCredit] = (1.0 - recovery_rates[iCredit]) / num_credits
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

    lossUnits = np.zeros(num_credits)
    numLossUnits = 1  # this is the zero loss

    for iCredit in range(0, num_credits):
        lossUnits[iCredit] = lossAmounts[iCredit] / gcd
        numLossUnits = numLossUnits + lossUnits[iCredit]

    defaultProbs = np.zeros(num_credits)

    for iCredit in range(0, num_credits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    lossDbn = lossDbnRecursionGCD(num_credits,
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
                       num_credits,
                       survivalProbabilities,
                       recovery_rates,
                       betaVector,
                       numIntegrationSteps):
    """ Get the approximated tranche survival probability of a portfolio
    of credits in the one-factor GC model using a Gaussian fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. Note that the losses in this fit are allowed to be negative. """

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    defaultProbs = [0.0] * num_credits
    for iCredit in range(0, num_credits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    dz = 2.0 * abs(minZ) / numIntegrationSteps
    z = minZ

    thresholds = np.zeros(num_credits)
    losses = np.zeros(num_credits)

    for iCredit in range(0, num_credits):
        pd = 1.0 - survivalProbabilities[iCredit]
        thresholds[iCredit] = norminvcdf(pd)
        losses[iCredit] = (1.0 - recovery_rates[iCredit]) / num_credits

    v = 0.0
    for _ in range(0, numIntegrationSteps):

        mu = 0.0
        var = 0.0

        # calculate the mean and variance of the conditional loss distribution
        for iCredit in range(0, num_credits):
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
def lossDbnHeterogeneousAdjBinomial(num_credits,
                                    defaultProbs,
                                    lossRatio,
                                    betaVector,
                                    numIntegrationSteps):
    """ Get the portfolio loss distribution using the adjusted binomial
    approximation to the conditional loss distribution. """

    numLossUnits = num_credits + 1
    condDefaultProbs = np.zeros(num_credits)
    uncondLossDbn = np.zeros(numLossUnits)
    indepDbn = np.zeros(numLossUnits)

    # Determine default threshold for each credit
    thresholds = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        thresholds[iCredit] = norminvcdf(defaultProbs[iCredit])

    dz = 2.0 * abs(minZ) / numIntegrationSteps
    z = minZ

    for iStep in range(0, numIntegrationSteps):
        for iCredit in range(0, num_credits):
            beta = betaVector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indepLossDbnHeterogeneousAdjBinomial(num_credits,
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
                          num_credits,
                          survivalProbabilities,
                          recovery_rates,
                          betaVector,
                          numIntegrationSteps):
    """ Get the approximated tranche survival probability of a portfolio of
    credits in the one-factor GC model using the adjusted binomial fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. This approach is both fast and highly accurate. """

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    defaultProbs = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        defaultProbs[iCredit] = 1.0 - survivalProbabilities[iCredit]

    totalLoss = 0.0
    for iCredit in range(0, num_credits):
        totalLoss += (1.0 - recovery_rates[iCredit])
    totalLoss /= num_credits

    avgLoss = totalLoss / num_credits

    lossRatio = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        lossRatio[iCredit] = (
            1.0 - recovery_rates[iCredit]) / num_credits / avgLoss

    lossDbn = lossDbnHeterogeneousAdjBinomial(num_credits,
                                              defaultProbs,
                                              lossRatio,
                                              betaVector,
                                              numIntegrationSteps)
    trancheEL = 0.0
    numLossUnits = num_credits + 1
    for iLossUnit in range(0, numLossUnits):
        loss = iLossUnit * avgLoss
        trancheLoss = min(loss, k2) - min(loss, k1)
        trancheEL += trancheLoss * lossDbn[iLossUnit]

    q = 1.0 - trancheEL / (k2 - k1)
    return q

###############################################################################
