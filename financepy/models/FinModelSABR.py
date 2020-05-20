##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import log, sqrt
from numba import njit, float64

##########################################################################


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        float64),
    fastmath=True,
    cache=True)
def blackVolFromSABR(alpha, beta, rho, nu, f, k, t):

    if abs(rho) >= 0.999999999:
        raise ValueError("Rho is a correlation and must be less than 1.0")

    b = 1.0 - beta
    fk = f * k
    m = f / k
    sigma = 0.1

    if abs(m - 1.0) > 1e-6:
        sigma = 1.0
        numTerm1 = ((alpha * b)**2.0) / (fk**b) / 24.0
        numTerm2 = rho * beta * nu * alpha / (fk**(b / 2.0)) / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        logM = log(m)
        z = nu / alpha * (fk**(b / 2.0)) * logM
        denom = (fk**(b / 2)) * (1.0 + (b**2) / 24.0 *
                                 (logM**2) + (b**4) / 1920.0 * (logM**4))
        x = log((sqrt(1.0 - 2.0 * rho * z + z * z) + z - rho) / (1.0 - rho))
        sigma = num * z / (denom * x)
    else:
        numTerm1 = ((alpha * b)**2) / (f**(2.0 * b)) / 24.0
        numTerm2 = rho * beta * nu * alpha / (f**b) / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        denom = f**b
        sigma = num / denom

    if sigma <= 0.0:
        raise ValueError("SABR Volatility <= 0%.")

    return sigma

##########################################################################
