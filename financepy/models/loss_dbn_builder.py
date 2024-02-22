##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from numba import njit, float64, int64
import numpy as np

from ..utils.math import pair_gcd

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_hetero_adj_binomial(num_credits,
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

    v_approx = 0.0
    v_exact = 0.0

    for i_credit in range(0, num_credits):
        loss_ratio2 = loss_ratio[i_credit] ** 2
        v_approx += loss_ratio2 * p * (1.0 - p)
        v_exact += loss_ratio2 * cond_probs[i_credit] \
            * (1.0 - cond_probs[i_credit])

    ###########################################################################

    mean_loss = p * num_credits
    mean_above = round(mean_loss + 1)
    mean_below = round(mean_loss)

    if mean_above > num_credits:
        mean_above = num_credits

    diff_above = mean_above - mean_loss
    diff_below = mean_below - mean_loss

    # DOK - TO DO - SIMPLIFY THIS CODE AS PER JOD PAPER
    term = diff_above * diff_above + \
        (diff_below * diff_below - diff_above * diff_above) * diff_above

    numer = v_exact - term
    denom = v_approx - term

    if abs(denom) < 1e-30:
        denom = 1e-30

    alpha = numer / denom
    epsilon_below = (1.0 - alpha) * diff_above
    epsilon_above = (1.0 - alpha) - epsilon_below

    for i_loss_unit in range(0, num_losses):
        indep_dbn[i_loss_unit] *= alpha

    indep_dbn[int(mean_below)] += epsilon_below
    indep_dbn[int(mean_above)] += epsilon_above
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

    portfolio_gcd = float(temp / scaling)
    return portfolio_gcd

###############################################################################


@njit(float64[:](int64, float64[:], float64[:]), fastmath=True, cache=True)
def indep_loss_dbn_recursion_gcd(num_credits,
                                 cond_default_probs,
                                 loss_units):

    num_loss_units = 1
    for i in range(0, len(loss_units)):
        num_loss_units += int(loss_units[i])

    prev_dbn = np.zeros(num_loss_units)
    prev_dbn[0] = 1.0

    small = 1e-10
    next_dbn = np.zeros(num_loss_units)

    for i_credit in range(0, num_credits):

        p = cond_default_probs[i_credit]
        loss = (int)(loss_units[i_credit] + small)

        for i_loss_unit in range(0, loss):
            next_dbn[i_loss_unit] = prev_dbn[i_loss_unit] * (1.0 - p)

        for i_loss_unit in range(loss, num_loss_units):
            next_dbn[i_loss_unit] = prev_dbn[i_loss_unit - loss] * \
                p + prev_dbn[i_loss_unit] * (1.0 - p)

        for i_loss_unit in range(0, num_loss_units):
            prev_dbn[i_loss_unit] = next_dbn[i_loss_unit]

    return next_dbn

##########################################################################
