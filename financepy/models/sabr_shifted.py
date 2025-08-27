# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import njit
from scipy.optimize import minimize

from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from ..utils.math import normcdf
from ..utils.helpers import label_to_string

# TODO: Should I merge this with SABR ?

########################################################################################


@njit(fastmath=True, cache=True)
def _x(rho, z):
    """Return function x used in Hagan's 2002 SABR lognormal vol expansion."""
    a = (1.0 - 2.0 * rho * z + z**2) ** 0.5 + z - rho
    b = 1.0 - rho
    return np.log(a / b)


########################################################################################


@njit(fastmath=True, cache=True)
def vol_function_shifted_sabr(params, f, k, t):
    """Black volatility implied by SABR model."""

    alpha = params[0]
    beta = params[1]
    rho = params[2]
    nu = params[3]
    shift = params[4]

    f = f + shift
    k = k + shift

    alpha = max(alpha, 1e-10)

    # Negative strikes or forwards
    if k <= 0:
        raise FinError("Strike must be positive")

    if f <= 0:
        raise FinError("Forward must be positive")

    logfk = np.log(f / k)
    b = 1.0 - beta
    fkb = (f * k) ** b
    a = b**2 * alpha**2 / (24.0 * fkb)
    b = 0.25 * rho * beta * nu * alpha / fkb**0.5
    c = (2.0 - 3.0 * rho**2.0) * nu**2.0 / 24
    d = fkb**0.5
    v = b**2 * logfk**2 / 24.0
    w = b**4 * logfk**4 / 1920.0
    z = nu * fkb**0.5 * logfk / alpha

    eps = 1e-07

    if abs(z) > eps:
        vz = alpha * z * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w) * _x(rho, z))
        return vz

    v0 = alpha * (1.0 + (a + b + c) * t) / (d * (1.0 + v + w))
    return v0


########################################################################################


class SABRShifted:
    """SABR - Shifted Stochastic alpha beta rho model by Hagan et al. is a
    stochastic volatility model where alpha controls the implied volatility,
    beta is the exponent on the the underlying asset's process so beta = 0
    is normal and beta = 1 is lognormal, rho is the correlation between the
    underlying and the volatility process. The shift allows negative rates."""

    ####################################################################################

    def __init__(self, alpha, beta, rho, nu, shift):
        """Create SABRShifted with all of the model parameters. We
        also provide functions below to assist with the calibration of the
        value of alpha."""

        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu
        self._shift = shift

    ####################################################################################

    def black_vol(self, f, k, t):
        """Black volatility from SABR model using Hagan et al. approx."""

        params = np.array([self._alpha, self._beta, self._rho, self._nu, self._shift])

        # I wish to enable vectorisations
        if isinstance(f, np.ndarray):
            vols = []
            for x in f:
                v = vol_function_shifted_sabr(params, x, k, t)
                vols.append(v)
            return np.array(vols)
        elif isinstance(k, np.ndarray):
            vols = []
            for x in k:
                v = vol_function_shifted_sabr(params, f, x, t)
                vols.append(v)
            return np.array(vols)
        elif isinstance(t, np.ndarray):
            vols = []
            for x in t:
                v = vol_function_shifted_sabr(params, f, k, x)
                vols.append(v)
            return np.array(vols)

        v = vol_function_shifted_sabr(params, f, k, t)
        return v

    ####################################################################################

    def black_vol_with_alpha(self, alpha, f, k, t):

        self._alpha = alpha[0]
        black_vol = self.black_vol(f, k, t)
        return black_vol

    ####################################################################################

    def value(
        self,
        forward_rate,  # Forward rate F
        strike_rate,  # Strike Rate k
        time_to_expiry,  # Time to Expiry (years)
        df,  # Discount Factor to expiry date
        call_or_put,
    ):  # Call or put
        """Price an option using Black's model which values in the forward
        measure following a change of measure."""

        f = forward_rate
        t = time_to_expiry
        k = strike_rate
        sqrt_t = np.sqrt(t)
        vol = self.black_vol(f, k, t)

        d1 = (np.log((f) / (k)) + vol * vol * t / 2) / (vol * sqrt_t)
        d2 = d1 - vol * sqrt_t

        if call_or_put == OptionTypes.EUROPEAN_CALL:
            return df * (f * normcdf(d1) - k * normcdf(d2))
        elif call_or_put == OptionTypes.EUROPEAN_PUT:
            return df * (k * normcdf(-d2) - f * normcdf(-d1))

        raise FinError("Option type must be a European Call(C) or Put(P)")

    ####################################################################################

    def set_alpha_from_black_vol(self, black_vol, forward, strike, time_to_expiry):
        """Estimate the valu normcdf(f the alpha coefficient of the SABR model
        by solving for the value of alpha that makes the SABR black vol equal
        to the input black vol. This uses a numerical 1D solver."""

        t_exp = time_to_expiry
        f = forward
        k = strike

        # The starting point is based on assuming that the strike is ATM
        self.set_alpha_from_atm_black_vol(black_vol, strike, time_to_expiry)

        init_alpha = self._alpha

        if init_alpha != black_vol:
            # Objective function

            ############################################################################

            def fn(x):

                return np.sqrt(
                    (black_vol - self.black_vol_with_alpha(x, f, k, t_exp)) ** 2
                )

            bnds = ((0.0, None),)
            x0 = init_alpha
            results = minimize(fn, x0, method="L-BFGS-B", bounds=bnds, tol=1e-8)
            alpha = results.x[0]
        else:
            alpha = init_alpha

        self._alpha = alpha

    ####################################################################################

    def set_alpha_from_atm_black_vol(self, black_vol, atm_strike, time_to_expiry):
        """We solve cubic equation for the unknown variable alpha for the
        special ATM case of the strike equalling the forward following Hagan
        and al. equation (3.3). We take the smallest real root as the preferred
        solution. This is useful for calibrating the model when beta has been
        chosen."""

        # For shifted SABR
        atm_strike = atm_strike + self._shift

        beta = self._beta
        rho = self._rho
        nu = self._nu
        t_exp = time_to_expiry
        k = atm_strike

        coeff0 = -black_vol * (k ** (1.0 - self._beta))
        coeff1 = 1.0 + ((2.0 - 3.0 * rho**2) / 24.0) * (nu**2) * t_exp
        coeff2 = (rho * beta * nu * t_exp) / (4.0 * (k ** (1.0 - beta)))
        coeff3 = (((1.0 - beta) ** 2) * t_exp) / (24.0 * (k ** (2.0 - 2.0 * beta)))
        coeffs = [coeff3, coeff2, coeff1, coeff0]
        roots = np.roots(coeffs)

        # Selecting the smallest positive real root
        alpha = np.min([coeff.real for coeff in roots if coeff.real > 0])
        self._alpha = alpha

    ####################################################################################

    def __repr__(self):
        """Return string with class details."""

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("Alpha", self._alpha)
        s += label_to_string("Beta", self._beta)
        s += label_to_string("Nu", self._nu)
        s += label_to_string("Rho", self._rho)
        s += label_to_string("Shift", self._shift)
        return s
