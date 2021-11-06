##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

##########################################################################

from ..utils.math import norminvcdf, N, INVROOT2PI
from ..utils.error import FinError
from .loss_dbn_builder import indep_loss_dbn_recursion_gcd
from .loss_dbn_builder import indep_loss_dbn_heterogeneous_adj_binomial
from .loss_dbn_builder import portfolio_gcd

###############################################################################

minZ = -6.0

###############################################################################
# This implements the one-factor latent variable formulation of the Gaussian
# Copula model as well as some approximations
###############################################################################


@njit(float64[:](int64, float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def loss_dbn_recursion_gcd(num_credits,
                           default_probs,
                           lossUnits,
                           beta_vector,
                           num_integration_steps):
    """ Full construction of the loss distribution of a portfolio of credits
    where losses have been calculate as number of units based on the GCD. """

    if len(default_probs) != num_credits:
        raise FinError("Default probability length must equal num credits.")

    if len(lossUnits) != num_credits:
        raise FinError("Loss units length must equal num credits.")

    if len(beta_vector) != num_credits:
        raise FinError("Beta vector length must equal num credits.")

    numLossUnits = 1
    for i in range(0, len(lossUnits)):
        numLossUnits += int(lossUnits[i])

    uncondLossDbn = np.zeros(numLossUnits)

    z = minZ
    dz = 2.0 * abs(z) / num_integration_steps

    condDefaultProbs = np.zeros(num_credits)

    thresholds = np.zeros(num_credits)
    for iCredit in range(0, int(num_credits)):
        thresholds[iCredit] = norminvcdf(default_probs[iCredit])

    for _ in range(0, num_integration_steps):
        for iCredit in range(0, num_credits):
            beta = beta_vector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indep_loss_dbn_recursion_gcd(num_credits,
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
def homog_basket_loss_dbn(survival_probabilities,
                          recovery_rates,
                          beta_vector,
                          num_integration_steps):
    """ Calculate the loss distribution of a CDS default basket where the
    portfolio is equally weighted and the losses in the portfolio are homo-
    geneous i.e. the credits have the same recovery rates. """

    num_credits = len(survival_probabilities)

    if num_credits == 0:
        raise FinError("Number of credits equals zero")

    for iCredit in range(1, num_credits):
        if recovery_rates[iCredit] != recovery_rates[0]:
            raise FinError("Losses are not homogeneous")

    m = 0.0
    for i in range(0, len(beta_vector)):
        m += beta_vector[i]
    m /= len(beta_vector)

    # High beta requires more integration steps
    if m > 0.7:
        num_integration_steps *= 2

    if m > 0.9:
        num_integration_steps *= 5

    lossUnits = np.ones(num_credits)

    default_probs = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        default_probs[iCredit] = 1.0 - survival_probabilities[iCredit]

    lossDbn = loss_dbn_recursion_gcd(num_credits,
                                     default_probs,
                                     lossUnits,
                                     beta_vector,
                                     num_integration_steps)

    return lossDbn

###############################################################################


@njit(float64(float64, float64, int64, float64[:], float64[:], float64[:],
              int64), fastmath=True)
def tranche_surv_prob_recursion(k1,
                                k2,
                                num_credits,
                                survival_probabilities,
                                recovery_rates,
                                beta_vector,
                                num_integration_steps):
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
    for i in range(0, len(beta_vector)):
        m += beta_vector[i]
    m /= len(beta_vector)

    if m > 0.8:
        num_integration_steps *= 2

    if commonRecoveryFlag == 1:
        gcd = lossAmounts[0]
    else:
        gcd = portfolio_gcd(lossAmounts)

    lossUnits = np.zeros(num_credits)
    numLossUnits = 1  # this is the zero loss

    for iCredit in range(0, num_credits):
        lossUnits[iCredit] = lossAmounts[iCredit] / gcd
        numLossUnits = numLossUnits + lossUnits[iCredit]

    default_probs = np.zeros(num_credits)

    for iCredit in range(0, num_credits):
        default_probs[iCredit] = 1.0 - survival_probabilities[iCredit]

    lossDbn = loss_dbn_recursion_gcd(num_credits,
                                     default_probs,
                                     lossUnits,
                                     beta_vector,
                                     num_integration_steps)

    trancheEL = 0.0
    for iLossUnit in range(0, int(numLossUnits)):
        loss = iLossUnit * gcd
        trancheLoss = min(loss, k2) - min(loss, k1)
        trancheEL = trancheEL + trancheLoss * lossDbn[iLossUnit]

    q = 1.0 - trancheEL / (k2 - k1)
    return q

###############################################################################


@njit(float64(float64, float64, float64, float64), fastmath=True, cache=True)
def gauss_approx_tranche_loss(k1, k2, mu, sigma):

    if abs(sigma) < 1e-6:
        gauss_approx_tranche_loss = 0.0
        if mu > k1:
            gauss_approx_tranche_loss += (mu - k1)

        if mu > k2:
            gauss_approx_tranche_loss += (mu - k2)
    else:
        d1 = (mu - k1) / sigma
        d2 = (mu - k2) / sigma

        gauss_approx_tranche_loss = (mu - k1) * N(d1) - (mu - k2) * N(d2)
        + sigma * np.exp(-0.5 * d1 * d1) * INVROOT2PI
        - sigma * np.exp(-0.5 * d2 * d2) * INVROOT2PI

    return gauss_approx_tranche_loss

###############################################################################


@njit(float64(float64, float64, int64, float64[:], float64[:], float64[:],
              int64), fastmath=True, cache=True)
def tranch_surv_prob_gaussian(k1,
                              k2,
                              num_credits,
                              survival_probabilities,
                              recovery_rates,
                              beta_vector,
                              num_integration_steps):
    """ Get the approximated tranche survival probability of a portfolio
    of credits in the one-factor GC model using a Gaussian fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. Note that the losses in this fit are allowed to be negative. """

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    default_probs = [0.0] * num_credits
    for iCredit in range(0, num_credits):
        default_probs[iCredit] = 1.0 - survival_probabilities[iCredit]

    dz = 2.0 * abs(minZ) / num_integration_steps
    z = minZ

    thresholds = np.zeros(num_credits)
    losses = np.zeros(num_credits)

    for iCredit in range(0, num_credits):
        pd = 1.0 - survival_probabilities[iCredit]
        thresholds[iCredit] = norminvcdf(pd)
        losses[iCredit] = (1.0 - recovery_rates[iCredit]) / num_credits

    v = 0.0
    for _ in range(0, num_integration_steps):

        mu = 0.0
        var = 0.0

        # calculate the mean and variance of the conditional loss distribution
        for iCredit in range(0, num_credits):
            beta = beta_vector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condprob = N(argz)
            mu += condprob * losses[iCredit]
            var += (losses[iCredit]**2) * condprob * (1.0 - condprob)

        sigma = np.sqrt(var)
        el = gauss_approx_tranche_loss(k1, k2, mu, sigma)
        gaussWt = np.exp(-(z**2) / 2.0)

        v += el * gaussWt
        z += dz

    v *= INVROOT2PI * dz
    q = 1.0 - v / (k2 - k1)
    return q

###############################################################################


@njit(float64[:](int64, float64[:], float64[:], float64[:], int64),
      fastmath=True, cache=True)
def loss_dbn_hetero_adj_binomial(num_credits,
                                 default_probs,
                                 loss_ratio,
                                 beta_vector,
                                 num_integration_steps):
    """ Get the portfolio loss distribution using the adjusted binomial
    approximation to the conditional loss distribution. """

    numLossUnits = num_credits + 1
    condDefaultProbs = np.zeros(num_credits)
    uncondLossDbn = np.zeros(numLossUnits)
    indepDbn = np.zeros(numLossUnits)

    # Determine default threshold for each credit
    thresholds = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        thresholds[iCredit] = norminvcdf(default_probs[iCredit])

    dz = 2.0 * abs(minZ) / num_integration_steps
    z = minZ

    for iStep in range(0, num_integration_steps):
        for iCredit in range(0, num_credits):
            beta = beta_vector[iCredit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[iCredit] - beta * z) / denom
            condDefaultProbs[iCredit] = N(argz)

        indepDbn = indep_loss_dbn_heterogeneous_adj_binomial(num_credits,
                                                             condDefaultProbs,
                                                             loss_ratio)

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
def tranche_surv_prob_adj_binomial(k1,
                                   k2,
                                   num_credits,
                                   survival_probabilities,
                                   recovery_rates,
                                   beta_vector,
                                   num_integration_steps):
    """ Get the approximated tranche survival probability of a portfolio of
    credits in the one-factor GC model using the adjusted binomial fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. This approach is both fast and highly accurate. """

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("K1 >= K2")

    default_probs = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        default_probs[iCredit] = 1.0 - survival_probabilities[iCredit]

    totalLoss = 0.0
    for iCredit in range(0, num_credits):
        totalLoss += (1.0 - recovery_rates[iCredit])
    totalLoss /= num_credits

    avgLoss = totalLoss / num_credits

    loss_ratio = np.zeros(num_credits)
    for iCredit in range(0, num_credits):
        loss_ratio[iCredit] = (
            1.0 - recovery_rates[iCredit]) / num_credits / avgLoss

    lossDbn = loss_dbn_hetero_adj_binomial(num_credits,
                                           default_probs,
                                           loss_ratio,
                                           beta_vector,
                                           num_integration_steps)
    trancheEL = 0.0
    numLossUnits = num_credits + 1
    for iLossUnit in range(0, numLossUnits):
        loss = iLossUnit * avgLoss
        trancheLoss = min(loss, k2) - min(loss, k1)
        trancheEL += trancheLoss * lossDbn[iLossUnit]

    q = 1.0 - trancheEL / (k2 - k1)
    return q

###############################################################################
