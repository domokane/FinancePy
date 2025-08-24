# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from numba import njit
import numpy as np


from ..utils.math import normcdf
from ..utils.math import normpdf, norminvcdf, M
from ..utils.error import FinError

########################################################################################


@njit(fastmath=True, cache=True)
def tr_surv_prob_lhp(k1, k2, num_credits, survival_probs, recovery_rates, beta):
    """Get the approximated tranche survival probability of a portfolio of
    credits in the one-factor GC model using the large portfolio limit which
    assumes a homogenous portfolio with an infinite number of credits. This
    approach is very fast but not so as accurate as the adjusted binomial."""

    if k1 == 0.0 and k2 == 0.0:
        return 0.0

    if k1 >= k2:
        raise FinError("k_1 >= k_2")

    p = 0.0
    portfolio_el = 0.0
    for i_credit in range(0, num_credits):
        pd = 1.0 - survival_probs[i_credit]
        p += pd
        portfolio_el += pd * (1.0 - recovery_rates[i_credit])

    if p == 0.0:
        return 1.0

    p = p / num_credits
    portfolio_el = portfolio_el / num_credits

    recovery = 1.0 - portfolio_el / p
    elk1 = exp_min_lk(k1, p, recovery, 1.0, beta)
    elk2 = exp_min_lk(k2, p, recovery, 1.0, beta)
    value = 1.0 - (elk2 - elk1) / (k2 - k1)
    return value


########################################################################################


@njit(fastmath=True, cache=True)
def portfolio_cdf_lhp(k, num_credits, qvector, recovery_rates, beta):

    p = 0.0
    portfolio_el = 0.0

    for j in range(0, num_credits):
        p += 1.0 - qvector[j]
        portfolio_el += (1.0 - recovery_rates[j]) * (1 - qvector[j])

    p = p / num_credits
    portfolio_el /= num_credits

    if p == 0:
        return 0.0

    recovery = 1.0 - portfolio_el / p

    if beta == 0:
        beta = 0.0000000001

    if k >= (1.0 - recovery):
        return 1.0

    if abs(beta) > 1.0:
        return 0.0

    c = norminvcdf(p)
    arga = k / (1.0 - recovery)
    a = 1.0 / beta * (c - np.sqrt(1.0 - beta * beta) * norminvcdf(arga))
    return normcdf(-a)


########################################################################################


@njit(fastmath=True, cache=True)
def exp_min_lk(k, p, r, n, beta):

    if beta == 0:
        beta = 0.0000000001

    if abs(beta) > 1.0:
        return 0.0

    if p == 0.0:
        return 0.0

    if k == 0.0:
        return 0.0

    if k >= (1.0 - r):
        return p * (1.0 - r)

    c = norminvcdf(p)
    arga = k / (1.0 - r) / n

    a = 1.0 / beta * (c - np.sqrt(1.0 - beta * beta) * norminvcdf(arga))
    el1 = (1.0 - r) * M(c, -a, -beta) + k * normcdf(a)
    return el1


########################################################################################


@njit(fastmath=True, cache=True)
def lhp_density(k, p, r, beta):

    if beta == 0.0:
        beta = 0.0000000001

    if abs(beta) > 0.99999:
        return 999

    if k == 0.0:
        return 0.0

    if k >= (1 - r):
        return 0.0

    c = norminvcdf(p)
    arga = k / (1.0 - r)

    dk = 0.0000001

    a = 1.0 / beta * (c - np.sqrt(1.0 - beta * beta) * norminvcdf(arga))
    term1 = normcdf(a)

    arga = (k + dk) / (1.0 - r)
    a = 1.0 / beta * (c - np.sqrt(1.0 - beta * beta) * norminvcdf(arga))

    term2 = normcdf(a)
    rho = -(term2 - term1) / dk

    return rho


########################################################################################


@njit(fastmath=True, cache=True)
def lhp_analytical_density_base_corr(k, p, r, beta, dbeta_dk):

    if beta == 0:
        beta = 0.0000000001

    if abs(beta) > 0.99999:
        return 99

    if k == 0.0:
        return 0.0

    if k >= (1.0 - r):
        return 0.0

    c = norminvcdf(p)
    arga = k / (1.0 - r)
    root_1_minus_beta_sqd = np.sqrt(1.0 - beta * beta)
    a = 1.0 / beta * (c - root_1_minus_beta_sqd * norminvcdf(arga))

    da_dk = -c / beta / beta * dbeta_dk
    da_dk = da_dk + (
        1.0 / root_1_minus_beta_sqd + root_1_minus_beta_sqd / beta / beta
    ) * dbeta_dk * norminvcdf(arga)
    da_dk = da_dk - np.sqrt(1.0 - beta * beta) / beta / normpdf(
        norminvcdf(k / (1.0 - r))
    ) / (1.0 - r)

    rho = -normpdf(a) * da_dk

    return rho


########################################################################################


@njit(fastmath=True, cache=True)
def lhp_analytical_density(k, p, r, beta):

    if beta == 0.0:
        beta = 1e-8

    if abs(beta) > 0.99999:
        raise Exception("Beta is greater than 0.99999")

    if k == 0:
        return 0.0

    if k >= (1.0 - r):
        return 0.0

    c = norminvcdf(p)
    arga = k / (1.0 - r)
    a = 1.0 / beta * (c - np.sqrt(1.0 - beta * beta) * norminvcdf(arga))
    da_dk = (
        -np.sqrt(1.0 - beta * beta)
        / beta
        / normpdf(norminvcdf(k / (1.0 - r)))
        / (1.0 - r)
    )
    rho = -normpdf(a) * da_dk

    return rho


"""

########################################################################################

@njit(fastmath=True, cache=True)
def exp_min_lk(k, p, r, n, beta):

    if abs(beta) > 1.0:
        raise FinError("Beta > 100%")

    if abs(beta) < 1e-10:
        beta = 1e-10

    if p == 0.0:
        return 0.0

    if k == 0.0:
        return 0.0

    if k >= (1.0 - r):
        return p * (1.0 - r)

    c = normpdf(p)
    arga = k / (1.0 - r) / n

    a = (1.0 / beta) * (c - np.sqrt(1.0 - beta * beta) * normpdf(arga))
    el = (1.0 - r) * M(c, -a, -beta) + k * normcdf(a)

    return el
"""

########################################################################################


@njit(fastmath=True, cache=True)
def prob_l_greater_than_k(k, p, r, beta):

    c = normpdf(p)
    arga = k / (1.0 - r)
    a = (1.0 / beta) * (c - np.sqrt(1.0 - beta * beta) * normpdf(arga))
    prob = 1.0 - normcdf(a)
    return prob
