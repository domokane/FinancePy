# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from numba import njit, float64, int64
import numpy as np


from ..utils.math import norminvcdf, normcdf, INV_ROOT_2_PI
from ..utils.error import FinError
from .loss_dbn_builder import indep_loss_dbn_recursion_gcd
from .loss_dbn_builder import indep_loss_dbn_hetero_adj_binomial
from .loss_dbn_builder import portfolio_gcd


MIN_Z = -6.0

# This implements the one-factor latent variable formulation of the Gaussian
# Copula model as well as some approximations

########################################################################################


@njit(
    float64[:](int64, float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
)
def loss_dbn_recursion_gcd(
    num_credits, default_probs, loss_units, beta_vector, num_integration_steps
):
    """Full construction of the loss distribution of a portfolio of credits
    where losses have been calculate as number of units based on the GCD."""

    if len(default_probs) != num_credits:
        raise FinError("Default probability length must equal num credits.")

    if len(loss_units) != num_credits:
        raise FinError("Loss units length must equal num credits.")

    if len(beta_vector) != num_credits:
        raise FinError("Beta vector length must equal num credits.")

    num_loss_units = 1
    for lu in loss_units:
        num_loss_units += int(lu)

    uncond_loss_dbn = np.zeros(num_loss_units)

    z = MIN_Z
    dz = 2.0 * abs(z) / num_integration_steps

    cond_default_probs = np.zeros(num_credits)

    thresholds = np.zeros(num_credits)
    for i_credit in range(0, int(num_credits)):
        thresholds[i_credit] = norminvcdf(default_probs[i_credit])

    for _ in range(0, num_integration_steps):
        for i_credit in range(0, num_credits):
            beta = beta_vector[i_credit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[i_credit] - beta * z) / denom
            cond_default_probs[i_credit] = normcdf(argz)

        indep_dbn = indep_loss_dbn_recursion_gcd(
            num_credits, cond_default_probs, loss_units
        )

        gauss_wt = np.exp(-(z * z) / 2.0)

        for i_unit in range(0, num_loss_units):
            uncond_loss_dbn[i_unit] += indep_dbn[i_unit] * gauss_wt

        z += dz

    for i_unit in range(0, int(num_loss_units)):
        uncond_loss_dbn[i_unit] *= INV_ROOT_2_PI * dz

    return uncond_loss_dbn


########################################################################################


@njit(
    float64[:](float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
)
def homog_basket_loss_dbn(
    survival_probs, recovery_rates, beta_vector, num_integration_steps
):
    """Calculate the loss distribution of a CDS default basket where the
    portfolio is equally weighted and the losses in the portfolio are homo-
    geneous i.e. the credits have the same recovery rates."""

    num_credits = len(survival_probs)

    if num_credits == 0:
        raise FinError("Number of credits equals zero")

    for i_credit in range(1, num_credits):
        if recovery_rates[i_credit] != recovery_rates[0]:
            raise FinError("Losses are not homogeneous")

    m = 0.0
    for beta in beta_vector:
        m += beta
    m /= len(beta_vector)

    # High beta requires more integration steps
    if m > 0.7:
        num_integration_steps *= 2

    if m > 0.9:
        num_integration_steps *= 5

    loss_units = np.ones(num_credits)

    default_probs = np.zeros(num_credits)
    for i_credit in range(0, num_credits):
        default_probs[i_credit] = 1.0 - survival_probs[i_credit]

    loss_dbn = loss_dbn_recursion_gcd(
        num_credits,
        default_probs,
        loss_units,
        beta_vector,
        num_integration_steps,
    )

    return loss_dbn


########################################################################################


@njit(
    float64(float64, float64, int64, float64[:], float64[:], float64[:], int64),
    fastmath=True,
)
def tranche_surv_prob_recursion(
    k1,
    k2,
    num_credits,
    survival_probs,
    recovery_rates,
    beta_vector,
    num_integration_steps,
):
    """Get the tranche survival probability of a portfolio of credits in the
    one-factor GC model using a full recursion calculation of the loss
    distribution and survival probabilities to some time horizon."""

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("k_1 >= k_2")

    common_recovery_flag = 1

    loss_amounts = np.zeros(num_credits)
    for i_credit in range(0, num_credits):
        loss_amounts[i_credit] = (1.0 - recovery_rates[i_credit]) / num_credits
        if loss_amounts[i_credit] != loss_amounts[0]:
            common_recovery_flag = 0

    gcd = 0.0

    m = 0.0
    for i in range(0, len(beta_vector)):
        m += beta_vector[i]
    m /= len(beta_vector)

    if m > 0.8:
        num_integration_steps *= 2

    if common_recovery_flag == 1:
        gcd = loss_amounts[0]
    else:
        gcd = portfolio_gcd(loss_amounts)

    loss_units = np.zeros(num_credits)
    num_loss_units = 1  # this is the zero loss

    for i_credit in range(0, num_credits):
        loss_units[i_credit] = loss_amounts[i_credit] / gcd
        num_loss_units = num_loss_units + loss_units[i_credit]

    default_probs = np.zeros(num_credits)

    for i_credit in range(0, num_credits):
        default_probs[i_credit] = 1.0 - survival_probs[i_credit]

    loss_dbn = loss_dbn_recursion_gcd(
        num_credits,
        default_probs,
        loss_units,
        beta_vector,
        num_integration_steps,
    )

    tranche_el = 0.0
    for i_loss_unit in range(0, int(num_loss_units)):
        loss = i_loss_unit * gcd
        tranche_loss = min(loss, k2) - min(loss, k1)
        tranche_el = tranche_el + tranche_loss * loss_dbn[i_loss_unit]

    q = 1.0 - tranche_el / (k2 - k1)
    return q


########################################################################################


@njit(float64(float64, float64, float64, float64), fastmath=True, cache=True)
def gauss_approx_tranche_loss(k1, k2, mu, sigma):

    if abs(sigma) < 1e-6:
        tranche_loss = 0.0
        if mu > k1:
            tranche_loss += mu - k1

        if mu > k2:
            tranche_loss += mu - k2
    else:

        d1 = (mu - k1) / sigma
        d2 = (mu - k2) / sigma

        expd1 = np.exp(-0.5 * d1 * d1)
        expd2 = np.exp(-0.5 * d2 * d2)

        tranche_loss = (
            (mu - k1) * normcdf(d1)
            - (mu - k2) * normcdf(d2)
            + sigma * expd1 * INV_ROOT_2_PI
            - sigma * expd2 * INV_ROOT_2_PI
        )

    return tranche_loss


########################################################################################


@njit(
    float64(float64, float64, int64, float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
)
def tranch_surv_prob_gaussian(
    k1,
    k2,
    num_credits,
    survival_probs,
    recovery_rates,
    beta_vector,
    num_integration_steps,
):
    """Get the approximated tranche survival probability of a portfolio
    of credits in the one-factor GC model using a Gaussian fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. Note that the losses in this fit are allowed to be negative."""

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("k_1 >= k_2")

    default_probs = [0.0] * num_credits
    for i_credit in range(0, num_credits):
        default_probs[i_credit] = 1.0 - survival_probs[i_credit]

    dz = 2.0 * abs(MIN_Z) / num_integration_steps
    z = MIN_Z

    thresholds = np.zeros(num_credits)
    losses = np.zeros(num_credits)

    for i_credit in range(0, num_credits):
        pd = 1.0 - survival_probs[i_credit]
        thresholds[i_credit] = norminvcdf(pd)
        losses[i_credit] = (1.0 - recovery_rates[i_credit]) / num_credits

    v = 0.0
    for _ in range(0, num_integration_steps):

        mu = 0.0
        var = 0.0

        # calculate the mean and variance of the conditional loss distribution
        for i_credit in range(0, num_credits):
            beta = beta_vector[i_credit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[i_credit] - beta * z) / denom
            condprob = normcdf(argz)
            mu += condprob * losses[i_credit]
            var += (losses[i_credit] ** 2) * condprob * (1.0 - condprob)

        sigma = np.sqrt(var)
        el = gauss_approx_tranche_loss(k1, k2, mu, sigma)
        gauss_wt = np.exp(-(z**2) / 2.0)

        v += el * gauss_wt
        z += dz

    v *= INV_ROOT_2_PI * dz
    q = 1.0 - v / (k2 - k1)
    return q


########################################################################################


@njit(
    float64[:](int64, float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
)
def loss_dbn_hetero_adj_binomial(
    num_credits, default_probs, loss_ratio, beta_vector, num_integration_steps
):
    """Get the portfolio loss distribution using the adjusted binomial
    approximation to the conditional loss distribution."""

    num_loss_units = num_credits + 1
    cond_default_probs = np.zeros(num_credits)
    uncond_loss_dbn = np.zeros(num_loss_units)
    indep_dbn = np.zeros(num_loss_units)

    # Determine default threshold for each credit
    thresholds = np.zeros(num_credits)
    for i_credit in range(0, num_credits):
        thresholds[i_credit] = norminvcdf(default_probs[i_credit])

    dz = 2.0 * abs(MIN_Z) / num_integration_steps
    z = MIN_Z

    for _ in range(0, num_integration_steps):

        for i_credit in range(0, num_credits):
            beta = beta_vector[i_credit]
            denom = np.sqrt(1.0 - beta * beta)
            argz = (thresholds[i_credit] - beta * z) / denom
            cond_default_probs[i_credit] = normcdf(argz)

        indep_dbn = indep_loss_dbn_hetero_adj_binomial(
            num_credits, cond_default_probs, loss_ratio
        )

        gauss_wt = np.exp(-z * z / 2.0)

        for i_loss_unit in range(0, num_loss_units):
            uncond_loss_dbn[i_loss_unit] += indep_dbn[i_loss_unit] * gauss_wt

        z = z + dz

    for i_loss_unit in range(0, num_loss_units):
        uncond_loss_dbn[i_loss_unit] *= INV_ROOT_2_PI * dz

    return uncond_loss_dbn


########################################################################################


@njit(
    float64(float64, float64, int64, float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
)
def tranche_surv_prob_adj_binomial(
    k1,
    k2,
    num_credits,
    survival_probs,
    recovery_rates,
    beta_vector,
    num_integration_steps,
    ########################################################################################
):
    """Get the approximated tranche survival probability of a portfolio of
    credits in the one-factor GC model using the adjusted binomial fit of the
    conditional loss distribution and survival probabilities to some time
    horizon. This approach is both fast and highly accurate."""

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("k_1 >= k_2")

    default_probs = np.zeros(num_credits)
    for i_credit in range(0, num_credits):
        default_probs[i_credit] = 1.0 - survival_probs[i_credit]

    total_loss = 0.0
    for i_credit in range(0, num_credits):
        total_loss += 1.0 - recovery_rates[i_credit]
    total_loss /= num_credits

    avg_loss = total_loss / num_credits

    loss_ratio = np.zeros(num_credits)
    for i_credit in range(0, num_credits):
        loss_ratio[i_credit] = (1.0 - recovery_rates[i_credit]) / num_credits / avg_loss

    loss_dbn = loss_dbn_hetero_adj_binomial(
        num_credits,
        default_probs,
        loss_ratio,
        beta_vector,
        num_integration_steps,
    )
    tranche_el = 0.0
    num_loss_units = num_credits + 1
    for i_loss_unit in range(0, num_loss_units):
        loss = i_loss_unit * avg_loss
        tranche_loss = min(loss, k2) - min(loss, k1)
        tranche_el += tranche_loss * loss_dbn[i_loss_unit]

    q = 1.0 - tranche_el / (k2 - k1)
    return q
