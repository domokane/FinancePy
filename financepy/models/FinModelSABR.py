##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from numba import njit, float64

from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinMath import N
from ..finutils.FinError import FinError

###############################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def volFunctionSABR(params, f, k, t):

    alpha = params[0]
    beta = params[1]
    rho = params[2]
    nu = params[3]

    if abs(rho) >= 0.999999999:
        raise FinError("Rho is a correlation and must be less than 1.0")

    b = 1.0 - beta
    fk = f * k
    m = f / k

    if abs(m - 1.0) > 1e-6:
        sigma = 1.0
        numTerm1 = ((alpha * b)**2.0) / (fk**b) / 24.0
        numTerm2 = rho * beta * nu * alpha / (fk**(b / 2.0)) / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        logM = np.log(m)
        z = nu / alpha * (fk**(b / 2.0)) * logM
        denom = (fk**(b / 2)) * (1.0 + (b**2) / 24.0 *
                                 (logM**2) + (b**4) / 1920.0 * (logM**4))

        x = np.log((np.sqrt(1.0 - 2.0*rho*z + z**2.0) + z - rho)/(1.0 - rho))
        sigma = num*z/(denom*x)

    else:
        # when the option is at the money
        numTerm1 = ((alpha * b)**2) / (f**(2.0 * b)) / 24.0
        numTerm2 = rho * beta * nu * alpha / (f**b) / 4.0
        numTerm3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (numTerm1 + numTerm2 + numTerm3) * t)
        denom = f**b
        sigma = num / denom

    if sigma <= 0.0:
        raise FinError("SABR Volatility <= 0%.")

    return sigma

###############################################################################


class FinModelSABR():
    ''' SABR - Stochastic alpha beta rho model by Hagan et al. '''

    def __init__(self, alpha, beta, rho, nu):
        ''' Create FinModelSABR with model parameters.'''
        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu

###############################################################################

    def blackVol(self, f, k, t):
        ''' Black volatility from SABR model using Hagan et al. approx. '''

        params = np.array([self._alpha, self._beta, self._rho, self._nu])

        # I wish to enable vectorisations
        if isinstance(f, np.ndarray):
            vols = []
            for x in f:
                v = volFunctionSABR(params, x, k, t)
                vols.append(v)
            return np.array(vols)

        elif isinstance(k, np.ndarray):
            vols = []
            for x in k:
                v = volFunctionSABR(params, f, x, t)
                vols.append(v)
            return np.array(vols)

        elif isinstance(t, np.ndarray):
            vols = []
            for x in t:
                v = volFunctionSABR(params, f, k, x)
                vols.append(v)
            return np.array(vols)
        else:
            v = volFunctionSABR(params, f, k, t)
            return v

###############################################################################

    def value(self,
              forwardRate,   # Forward rate
              strikeRate,    # Strike Rate
              timeToExpiry,  # time to expiry in years
              df,            # Discount Factor to expiry date
              callOrPut):    # Call or put
        ''' Price an option using Black's model which values in the forward
        measure following a change of measure. '''

        f = forwardRate
        t = timeToExpiry
        k = strikeRate
        sqrtT = np.sqrt(t)
        vol = self.blackVol(f, k, t)

        d1 = np.log(f/k) + vol * vol * t / 2
        d1 = d1 / (vol * sqrtT)
        d2 = d1 - vol * sqrtT

        if callOrPut == FinOptionTypes.EUROPEAN_CALL:
            return df * (f * N(d1) - k * N(d2))
        elif callOrPut == FinOptionTypes.EUROPEAN_PUT:
            return df * (k * N(-d2) - f * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

        return 999

###############################################################################
