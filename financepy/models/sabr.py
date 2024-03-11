##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Volatility calibration functionality - Guillaume Lefieux 2021
###############################################################################

import numpy as np
from numba import njit, float64
from scipy.optimize import minimize

from ..utils.global_types import OptionTypes
from ..utils.math import N
from ..utils.error import FinError
from ..utils.helpers import label_to_string

###############################################################################


@njit
def _x(rho, z):
    """ Return function x used in Hagan's 2002 SABR lognormal vol expansion."""
    a = (1.0 - 2.0*rho*z + z**2)**.5 + z - rho
    b = 1.0 - rho
    return np.log(a / b)

##############################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def vol_function_sabr(params, f, k, t):
    """ Black volatility implied by SABR model. """

    alpha = params[0]
    beta = params[1]
    rho = params[2]
    nu = params[3]

    alpha = max(1e-10, alpha)

    # Negative strikes or forwards
    if k <= 0:
        raise FinError("Strike must be positive")

    if f <= 0:
        raise FinError("Forward must be positive")

    ln_f_over_k = np.log(f / k)
    b = 1.0 - beta
    fkb = (f*k)**b
    a = b**2 * alpha**2 / (24.0 * fkb)
    b = 0.25 * rho * beta * nu * alpha / fkb**0.5
    c = (2.0 - 3.0*rho**2.0) * nu**2.0 / 24
    d = fkb**0.5
    v = b**2 * ln_f_over_k**2 / 24.0
    w = b**4 * ln_f_over_k**4 / 1920.0
    z = nu * fkb**0.5 * ln_f_over_k / alpha

    eps = 1e-07

    if abs(z) > eps:
        vz = alpha * z * (1.0 + (a + b + c) * t) / \
            (d * (1.0 + v + w) * _x(rho, z))
        return vz
    else:
        v0 = alpha * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w))
        return v0

###############################################################################


@njit(float64(float64[:], float64, float64, float64),
      fastmath=True, cache=True)
def vol_function_sabr_beta_one(params, f, k, t):
    """ This is the SABR function with the exponent beta set equal to 1 so only
    3 parameters are free. The first parameter is alpha, then nu and the third
    parameter is rho. Check the order as it is not the same as main SABR fn"""

    alpha = params[0]
    rho = params[1]
    nu = params[2]

    if rho > 1.0:
        rho = 0.99

    if rho < -1.0:
        rho = -0.99

    m = f / k

    if abs(m - 1.0) > 1e-6:

        sigma = 1.0
        num_term1 = 0.0
        num_term2 = rho * nu * alpha / 4.0
        num_term3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (num_term1 + num_term2 + num_term3) * t)
        log_m = np.log(m)
        z = nu / alpha * log_m
        denom = 1.0
        x = np.log((np.sqrt(1.0 - 2.0*rho*z + z**2.0) + z - rho)/(1.0 - rho))
        sigma = num*z/(denom*x)

    else: # when the option is at the money

        num_term1 = 0.0
        num_term2 = rho * nu * alpha / 4.0
        num_term3 = nu * nu * ((2.0 - 3.0 * (rho**2.0)) / 24.0)
        num = alpha * (1.0 + (num_term1 + num_term2 + num_term3) * t)
        denom = 1.0
        sigma = num / denom

    return sigma

###############################################################################


@njit(float64(float64[:], float64, float64, float64),
      fastmath=True, cache=True)
def vol_function_sabr_beta_half(params, f, k, t):
    """ Black volatility implied by SABR model. """

    alpha = params[0]
    rho = params[1]
    nu = params[2]
    beta = 0.50

    alpha = max(1e-10, alpha)

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
        vz = alpha * z * (1.0 + (a + b + c) * t) / \
            (d * (1.0 + v + w) * _x(rho, z))
        return vz
    else:
        v0 = alpha * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w))
        return v0

###############################################################################


class SABR():
    """ SABR - Stochastic alpha beta rho model by Hagan et al. which is a
    stochastic volatility model where alpha controls the implied volatility,
    beta is the exponent on the the underlying asset's process so beta = 0
    is normal and beta = 1 is lognormal, rho is the correlation between the
    underlying and the volatility process. """

    def __init__(self, alpha, beta, rho, nu):
        """ Create SABR with all of the model parameters. We will
        also provide functions below to assist with the calibration of the
        value of alpha. """

        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.nu = nu

###############################################################################

    def black_vol(self, f, k, t):
        """ Black volatility from SABR model using Hagan et al. approx. """

        params = np.array([self.alpha, self.beta, self.rho, self.nu])

        # I wish to enable vectorisations
        if isinstance(f, np.ndarray):
            vols = []
            for x in f:
                v = vol_function_sabr(params, x, k, t)
                vols.append(v)
            return np.array(vols)

        elif isinstance(k, np.ndarray):
            vols = []
            for x in k:
                v = vol_function_sabr(params, f, x, t)
                vols.append(v)
            return np.array(vols)

        elif isinstance(t, np.ndarray):
            vols = []
            for x in t:
                v = vol_function_sabr(params, f, k, x)
                vols.append(v)
            return np.array(vols)
        
        v = vol_function_sabr(params, f, k, t)
        return v

###############################################################################

    def black_vol_with_alpha(self, alpha, f, k, t):

        self.alpha = alpha[0]
        black_vol = self.black_vol(f, k, t)
        return black_vol

###############################################################################

    def value(self,
              forward_rate,   # Forward rate
              strike_rate,    # Strike Rate
              time_to_expiry,  # time to expiry in years
              df,            # Discount Factor to expiry date
              call_or_put):    # Call or put
        """ Price an option using Black's model which values in the forward
        measure following a change of measure. """

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        sqrt_t = np.sqrt(t)
        vol = self.black_vol(f, k, t)

        d1 = np.log(f/k) + vol * vol * t / 2
        d1 = d1 / (vol * sqrt_t)
        d2 = d1 - vol * sqrt_t

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * (f * N(d1) - k * N(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * (k * N(-d2) - f * N(-d1))
        else:
            raise Exception("Option type must be a European Call(C) or Put(P)")

###############################################################################

    def set_alpha_from_black_vol(self,
                                 black_vol,
                                 forward,
                                 strike,
                                 time_to_expiry):
        """ Estimate the value of the alpha coefficient of the SABR model
        by solving for the value of alpha that makes the SABR black vol equal
        to the input black vol. This uses a numerical 1D solver. """
        t_exp = time_to_expiry
        f = forward
        k = strike

        # The starting point is based on assuming that the strike is ATM
        self.set_alpha_from_atm_black_vol(black_vol, strike, time_to_expiry)

        init_alpha = self.alpha

        if init_alpha != black_vol:
            # Objective function
            def fn(x): return np.sqrt(
                (black_vol - self.black_vol_with_alpha(x, f, k, t_exp)) ** 2)
            bnds = ((0.0, None),)
            x0 = init_alpha
            results = minimize(fn, x0, method="L-BFGS-B",
                               bounds=bnds, tol=1e-8)
            alpha = results.x[0]
        else:
            alpha = init_alpha

        self.alpha = alpha

###############################################################################

    def set_alpha_from_atm_black_vol(self, black_vol, atm_strike, time_to_expiry):
        """ We solve cubic equation for the unknown variable alpha for the
        special ATM case of the strike equalling the forward following Hagan
        and al. equation (3.3). We take the smallest real root as the preferred
        solution. This is useful for calibrating the model when beta has been
        chosen."""
        beta = self.beta
        rho = self.rho
        nu = self.nu
        t_exp = time_to_expiry
        K = atm_strike

        coeff0 = -black_vol * (K**(1.0 - self.beta))
        coeff1 = 1.0 + ((2.0 - 3.0 * rho**2) / 24.0) * (nu**2) * t_exp
        coeff2 = (rho * beta * nu * t_exp) / (4.0 * (K**(1.0 - beta)))
        coeff3 = (((1.0 - beta)**2) * t_exp) / (24.0 * (K**(2.0 - 2.0 * beta)))
        coeffs = [coeff3, coeff2, coeff1, coeff0]
        roots = np.roots(coeffs)

        # Selecting the smallest positive real root
        alpha = np.min([coeff.real for coeff in roots if coeff.real > 0])
        self.alpha = alpha

###############################################################################

    def __repr__(self):
        """ Return string with class details. """
        s = label_to_string("OBJECT TYPE", type(self)._name__)
        s += label_to_string("Alpha", self.alpha)
        s += label_to_string("Beta", self.beta)
        s += label_to_string("Nu", self.nu)
        s += label_to_string("Rho", self.rho)
        return s

###############################################################################
