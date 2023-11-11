##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

from ..utils.math import pair_gcd

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_heterogeneous_adj_binomial(num_credits,
                                              cond_probs,
                                              loss_ratio):

    # Algorithm due to D. O'Kane.

    num_losses = num_credits + 1
    indep_dbn = np.zeros(num_losses)

    p = 0.0
    for i_credit in range(0, num_credits):
        p += loss_ratio[i_credit] * cond_probs[i_credit]
    p = p / num_credits

    ###########################################################################

    if p < 0.5:
        ratio = p / (1.0 - p)
        indep_dbn[0] = (1.0 - p)**num_credits
        for i in range(1, num_losses):
            indep_dbn[i] = indep_dbn[i - 1] * ratio \
                * (num_credits - i + 1.0) / i
    else:
        ratio = (1.0 - p) / p
        indep_dbn[num_credits] = p ** num_credits
        for i in range(num_credits - 1, -1, -1):
            indep_dbn[i] = indep_dbn[i + 1] * \
                ratio * (i + 1.0) / (num_credits - i)

    ###########################################################################

    vapprox = 0.0
    vexact = 0.0

    for i_credit in range(0, num_credits):
        loss_ratio2 = loss_ratio[i_credit] ** 2
        vapprox += loss_ratio2 * p * (1.0 - p)
        vexact += loss_ratio2 * cond_probs[i_credit] \
            * (1.0 - cond_probs[i_credit])

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

    for i_loss_unit in range(0, num_losses):
        indep_dbn[i_loss_unit] *= alpha

    indep_dbn[int(meanBelow)] += epsilonBelow
    indep_dbn[int(meanAbove)] += epsilonAbove

    return indep_dbn

###############################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def portfolio_gcd(actual_losses):

    num_credits = len(actual_losses)
    scaling = 1000000

    temp = (int)(actual_losses[0] * scaling)

    for i_credit in range(1, num_credits):
        num2 = int(actual_losses[i_credit] * scaling)
        temp = pair_gcd(temp, num2)

    portfolioGCD = float(temp / scaling)
    return portfolioGCD

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_recursion_gcd(num_credits,
                                 cond_default_probs,
                                 loss_units):

    num_loss_units = 1
    for i in range(0, len(loss_units)):
        num_loss_units += int(loss_units[i])

    prevDbn = np.zeros(num_loss_units)
    prevDbn[0] = 1.0

    small = 1e-10
    nextDbn = np.zeros(num_loss_units)

    for i_credit in range(0, num_credits):

        p = cond_default_probs[i_credit]
        loss = (int)(loss_units[i_credit] + small)

        for i_loss_unit in range(0, loss):
            nextDbn[i_loss_unit] = prevDbn[i_loss_unit] * (1.0 - p)

        for i_loss_unit in range(loss, num_loss_units):
            nextDbn[i_loss_unit] = prevDbn[i_loss_unit - loss] * \
                p + prevDbn[i_loss_unit] * (1.0 - p)

        for i_loss_unit in range(0, num_loss_units):
            prevDbn[i_loss_unit] = nextDbn[i_loss_unit]

    return nextDbn

##########################################################################
