##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit, float64
from scipy.optimize import minimize

from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinError import FinError
from ..finutils.FinMath import N
from ..finutils.FinHelperFunctions import labelToString

###############################################################################
# TODO: Should I merge this with SABR ?
###############################################################################

###############################################################################
###############################################################################

@njit
def _x(rho, z):
    '''Return function x used in Hagan's 2002 SABR lognormal vol expansion.'''
    a = (1.0 - 2.0*rho*z + z**2)**.5 + z - rho
    b = 1.0 - rho
    return np.log(a / b)

@njit
def volFunctionShiftedSABR(params, f, k, t):
    ''' Black volatility implied by SABR model. '''

    alpha = params[0]
    beta = params[1]
    rho = params[2]
    nu = params[3]
    shift = params[4]

    f = f + shift
    k = k + shift

    if alpha < 1e-10:
        alpha = 1e-10

    # Negative strikes or forwards
    if k <= 0:
        raise FinError("Strike must be positive")

    if f <= 0:
        raise FinError("Forward must be positive")

    logfk = np.log(f / k)
    b = 1.0 - beta
    fkb = (f*k)**b
    a = b**2 * alpha**2 / (24.0 * fkb)
    b = 0.25 * rho * beta * nu * alpha / fkb**0.5
    c = (2.0 - 3.0*rho**2.0) * nu**2.0 / 24
    d = fkb**0.5
    v = b**2 * logfk**2 / 24.0
    w = b**4 * logfk**4 / 1920.0
    z = nu * fkb**0.5 * logfk / alpha

    eps = 1e-07

    if abs(z) > eps:
        vz = alpha * z * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w) * _x(rho, z))
        return vz
    else:
        v0 = alpha * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w))
        return v0

###############################################################################


class FinModelSABRShifted():
    ''' SABR - Shifted Stochastic alpha beta rho model by Hagan et al. is a 
    stochastic volatility model where alpha controls the implied volatility,
    beta is the exponent on the the underlying asset's process so beta = 0
    is normal and beta = 1 is lognormal, rho is the correlation between the 
    underlying and the volatility process. The shift allows negative rates.'''

    def __init__(self, alpha, beta, rho, nu, shift):
        ''' Create FinModelSABRShifted with all of the model parameters. We 
        also provide functions below to assist with the calibration of the 
        value of alpha. '''

        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu
        self._shift = shift

###############################################################################

    def blackVol(self, f, k, t):
        ''' Black volatility from SABR model using Hagan et al. approx. '''

        params = np.array([self._alpha, self._beta, self._rho, 
                           self._nu, self._shift])

        # I wish to enable vectorisations
        if isinstance(f, np.ndarray):
            vols = []
            for x in f:
                v = volFunctionShiftedSABR(params, x, k, t)
                vols.append(v)
            return np.array(vols)
        elif isinstance(k, np.ndarray):
            vols = []
            for x in k:
                v = volFunctionShiftedSABR(params, f, x, t)
                vols.append(v)
            return np.array(vols)

        elif isinstance(t, np.ndarray):
            vols = []
            for x in t:
                v = volFunctionShiftedSABR(params, f, k, x)
                vols.append(v)
            return np.array(vols)
        else:
            v = volFunctionShiftedSABR(params, f, k, t)
            return v

###############################################################################

    def blackVolWithAlpha(self, alpha, f, k, t):

        self._alpha = alpha[0]
        blackVol = self.blackVol(f, k, t)
        return blackVol

###############################################################################

    def value(self,
              forwardRate,   # Forward rate F
              strikeRate,    # Strike Rate K
              timeToExpiry,  # Time to Expiry (years)
              df,            # Discount Factor to expiry date
              callOrPut):    # Call or put
        ''' Price an option using Black's model which values in the forward
        measure following a change of measure. '''

        f = forwardRate
        t = timeToExpiry
        k = strikeRate
        sqrtT = np.sqrt(t)
        vol = self.blackVol(f, k, t)

        d1 = (np.log((f)/(k)) + vol * vol * t / 2) / (vol * sqrtT)
        d2 = d1 - vol*sqrtT

        if callOrPut == FinOptionTypes.EUROPEAN_CALL:
            return df * (f * N(d1) - k * N(d2))
        elif callOrPut == FinOptionTypes.EUROPEAN_PUT:
            return df * (k * N(-d2) - f * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

###############################################################################

    def setAlphaFromBlackVol(self, blackVol, forward, strike, timeToExpiry):
        ''' Estimate the value of the alpha coefficient of the SABR model
        by solving for the value of alpha that makes the SABR black vol equal
        to the input black vol. This uses a numerical 1D solver. '''

        texp = timeToExpiry
        f = forward
        K = strike

        # The starting point is based on assuming that the strike is ATM
        self.setAlphaFromATMBlackVol(blackVol, strike, timeToExpiry)

        initAlpha = self._alpha

        if initAlpha != blackVol:
            # Objective function
            fn = lambda x: np.sqrt((blackVol - self.blackVolWithAlpha(x, f, K, texp))**2)
            bnds = ((0.0, None),)
            x0 = initAlpha
            results = minimize(fn, x0, method="L-BFGS-B", bounds=bnds, tol=1e-8)
            alpha = results.x[0]
        else:
            alpha = initAlpha
        
        self._alpha = alpha

###############################################################################

    def setAlphaFromATMBlackVol(self, blackVol, atmStrike, timeToExpiry):
        ''' We solve cubic equation for the unknown variable alpha for the 
        special ATM case of the strike equalling the forward following Hagan 
        and al. equation (3.3). We take the smallest real root as the preferred
        solution. This is useful for calibrating the model when beta has been
        chosen.''' 

        # For shifted SABR
        atmStrike = atmStrike + self._shift

        beta = self._beta
        rho = self._rho
        nu = self._nu
        texp = timeToExpiry
        K = atmStrike

        coeff0 = -blackVol * (K**(1.0 - self._beta))
        coeff1 = 1.0 + ((2.0 - 3.0 * rho**2) / 24.0) * (nu**2) * texp
        coeff2 = (rho * beta * nu * texp) / (4.0 * (K**(1.0 - beta)))
        coeff3 = (((1.0 - beta)**2) * texp) / (24.0 * (K**(2.0 - 2.0 * beta)))
        coeffs = [coeff3, coeff2, coeff1, coeff0]
        roots = np.roots(coeffs)

        # Selecting the smallest positive real root
        alpha = np.min([coeff.real for coeff in roots if coeff.real > 0])
        self._alpha = alpha

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("Alpha", self._alpha)
        s += labelToString("Beta", self._beta)
        s += labelToString("Nu", self._nu)
        s += labelToString("Rho", self._rho)
        s += labelToString("Shift", self._shift)
        return s

###############################################################################
