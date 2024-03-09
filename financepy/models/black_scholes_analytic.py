##############################################################################
# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi
##############################################################################

import numpy as np
from numba import float64, int64, vectorize, njit

from ..utils.global_types import OptionTypes
from ..utils.global_vars import g_small
from ..utils.math import N, n_vect, n_prime_vect
from ..utils.error import FinError
from ..utils.solver_1d import bisection, newton, newton_secant

###############################################################################
# Analytical Black Scholes model implementation and approximations
###############################################################################


@vectorize([float64(float64, float64, float64, float64, float64, float64,
                    int64)], fastmath=True, cache=True)
def bs_value(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif option_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t

    value = phi * ss * N(phi * d1) - phi * kk * N(phi * d2)
    return value

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_delta(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""
    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = +1.0
    elif option_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Error: Unknown option type value")

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    delta = phi * np.exp(-q*t) * n_vect(phi * d1)
    return delta

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_gamma(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    gamma = np.exp(-q*t) * n_prime_vect(d1) / s / v_sqrt_t
    return gamma

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_vega(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""
    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    vega = ss * sqrt_t * n_prime_vect(d1)
    return vega

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_theta(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif option_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    theta = - ss * n_prime_vect(d1) * v / 2.0 / sqrt_t
    theta = theta - phi * r * k * np.exp(-r*t) * n_vect(phi * d2)
    theta = theta + phi * q * ss * n_vect(phi * d1)
    return theta

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_rho(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif option_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    rho = phi * k * t * np.exp(-r*t) * n_vect(phi * d2)
    return rho

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, float64, int64)], fastmath=True, cache=True)
def bs_vanna(s, t, k, r, q, v, option_type_value):
    """Price a derivative using Black-Scholes model."""

    k = np.maximum(k, g_small)
    t = np.maximum(t, g_small)
    v = np.maximum(v, g_small)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q*t)
    kk = k * np.exp(-r*t)
    d1 = np.log(ss/kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    vanna = np.exp(-q*t) * sqrt_t * n_prime_vect(d1) * (d2/v)
    return vanna

###############################################################################


@njit(fastmath=True, cache=True)
def _f(sigma, args):

    s = args[0]
    t = args[1]
    k = args[2]
    r = args[3]
    q = args[4]
    price = args[5]
    option_type_value = int(args[6])

    bs_price = bs_value(s, t, k, r, q, sigma, option_type_value)
    obj = bs_price - price
    return obj

##############################################################################


@njit(fastmath=True, cache=True)
def _fvega(sigma, args):

    s = args[0]
    t = args[1]
    k = args[2]
    r = args[3]
    q = args[4]
    option_type_value = int(args[6])
    vega = bs_vega(s, t, k, r, q, sigma, option_type_value)
    return vega

###############################################################################


@vectorize([float64(float64, float64, float64, float64,
                    float64, int64)], fastmath=True, cache=True)
def bs_intrinsic(s, t, k, r, q, option_type_value):
    """Calculate the Black-Scholes implied volatility of a European
     vanilla option using Newton with a fallback to bisection."""

    fwd = s * np.exp((r-q)*t)

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r*t) * max(fwd - k, 0.0)
    else:
        intrinsic_value = np.exp(-r*t) * max(k - fwd, 0.0)

    return intrinsic_value

###############################################################################


@vectorize([float64(float64, float64, float64, float64, float64, float64,
                    int64)], fastmath=True, cache=True,  forceobj=True)
def bs_implied_volatility(s, t, k, r, q, price, option_type_value):
    """Calculate the Black-Scholes implied volatility of a European
     vanilla option using Newton with a fallback to bisection."""

    fwd = s * np.exp((r-q)*t)

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r*t) * max(fwd - k, 0.0)
    else:
        intrinsic_value = np.exp(-r*t) * max(k - fwd, 0.0)

    div_adj_stock_price = s * np.exp(-q * t)
    df = np.exp(-r * t)

    # Flip ITM call option to be OTM put and vice-versa using put call parity
    if intrinsic_value > 0.0:

        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            price = price - (div_adj_stock_price - k * df)
            option_type_value = OptionTypes.EUROPEAN_PUT.value
        else:
            price = price + (div_adj_stock_price - k * df)
            option_type_value = OptionTypes.EUROPEAN_CALL.value

        # Update intrinsic based on new option type
        if option_type_value == OptionTypes.EUROPEAN_CALL.value:
            intrinsic_value = np.exp(-r*t) * max(fwd - k, 0.0)
        else:
            intrinsic_value = np.exp(-r*t) * max(k - fwd, 0.0)

    time_value = price - intrinsic_value

    # Add a tolerance in case it is just numerical imprecision
    if time_value < 0.0:
        print("Time value", time_value)
        raise FinError("Option Price is below the intrinsic value")

    ###########################################################################
    # Some approximations which might be used later
    ###########################################################################

    if option_type_value == OptionTypes.EUROPEAN_CALL.value:
        call = price
    else:
        call = price + (div_adj_stock_price - k * df)

    # Notation in SSRN-id567721.pdf
    X = k * np.exp(-r*t)
    S = s*np.exp(-q*t)
    pi = np.pi

    ###########################################################################
    # Initial point of inflexion
    ###########################################################################

    # arg = np.abs(np.log(fwd/k))
    # sigma0 = np.sqrt(2.0 * arg)

    ###########################################################################
    # Corrado Miller from Hallerbach equation (7)
    ###########################################################################

    cmsigma = 0.0
    # arg = (C - 0.5*(S-X))**2 - ((S-X)**2)/ pi

    # if arg < 0.0:
    #     arg = 0.0

    # cmsigma = (C-0.5*(S-X) + np.sqrt(arg))
    # cmsigma = cmsigma * np.sqrt(2.0*pi) / (S+X)
    # cmsigma = cmsigma / np.sqrt(t)

    ###########################################################################
    # Hallerbach SSRN-id567721.pdf equation (22)
    ###########################################################################

    hsigma = 0.0
    gamma = 2.0
    arg = (2*call+X-S)**2 - gamma * (S+X)*(S-X)*(S-X) / pi / S
    arg = max(arg, 0.0)

    hsigma = 2.0 * call + X - S + np.sqrt(arg)
    hsigma = hsigma * np.sqrt(2.0*pi) / 2.0 / (S+X)
    hsigma = hsigma / np.sqrt(t)

    sigma0 = hsigma

    ###########################################################################

    arglist = [s, t, k, r, q, price, option_type_value]
    argsv = np.array(arglist)

    tol = 1e-6
    sigma = newton(_f, sigma0, _fvega, argsv, tol=tol)

    if sigma is None:
        sigma = bisection(_f, 1e-4, 10.0, argsv, xtol=tol)
        if sigma is None:
            method = "Failed"
        else:
            method = "Bisection"
    else:
        method = "Newton"

    if 1 == 0:
        print("S: %7.2f K: %7.3f T:%5.3f V:%10.7f Sig0: %7.5f CM: %7.5f HL: %7.5f NW: %7.5f %10s" % (
            s, k, t, price, sigma0*100.0, cmsigma*100.0, hsigma*100.0, sigma*100.0, method))

    return sigma

###############################################################################
###############################################################################
# This module contains a number of analytical approximations for the price of
# an American style option starting with Barone-Adesi-Whaley
# https://deriscope.com/docs/Barone_Adesi_Whaley_1987.pdf
###############################################################################
###############################################################################


@njit(fastmath=True, cache=True)
def _fcall(si, *args):
    """Function to determine ststar for pricing American call options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    mm = 2.0 * r / v2
    ww = 2.0 * b / v2
    kk = 1.0 - np.exp(-r * t)

    q2 = (1.0 - ww + np.sqrt((ww - 1.0)**2 + 4.0 * mm/kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    obj_fn = si - k
    obj_fn = obj_fn - bs_value(si, t, k, r, q, v, +1)
    obj_fn = obj_fn - (1.0 - np.exp(-q*t) * n_vect(d1)) * si / q2
    return obj_fn

###############################################################################


@njit(fastmath=True, cache=True)
def _fput(si, *args):
    """Function to determine sstar for pricing American put options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    W = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q1 = (1.0 - W - np.sqrt((W - 1.0)**2 + 4.0 * K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    obj_fn = si - k
    obj_fn = obj_fn - bs_value(si, t, k, r, q, v, -1)
    obj_fn = obj_fn - (1.0 - np.exp(-q*t) * n_vect(-d1)) * si / q1
    return obj_fn

###############################################################################
# TODO: NUMBA SPEED UP
###############################################################################


@njit(fastmath=True)
def baw_value(s, t, k, r, q, v, phi):
    """American Option Pricing Approximation using the Barone-Adesi-Whaley
     approximation for the Black Scholes Model"""

    b = r - q

    if phi == 1:

        if b >= r:
            return bs_value(s, t, k, r, q, v, +1)

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fcall, x0=s, args=argtuple,
                              tol=1e-7, maxiter=50)

        M = 2.0 * r / (v*v)
        W = 2.0 * b / (v*v)
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar/k) + (b + v*v / 2.0) * t) / (v * np.sqrt(t))
        q2 = (-1.0 * (W - 1.0) + np.sqrt((W - 1.0)**2 + 4.0 * M/K)) / 2.0
        A2 = (sstar / q2) * (1.0 - np.exp(-q * t) * n_vect(d1))

        if s < sstar:
            return bs_value(s, t, k, r, q, v, OptionTypes.EUROPEAN_CALL.value) + A2 * ((s/sstar)**q2)

        return s - k

    elif phi == -1:

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fput, x0=s, args=argtuple, tol=1e-7, maxiter=50)

        v2 = v * v

        M = 2.0 * r / v2
        W = 2.0 * b / v2
        K = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (-1.0 * (W - 1.0) - np.sqrt((W - 1.0)**2 + 4.0 * M/K)) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * n_vect(-d1))

        if s > sstar:
            bsv = bs_value(s, t, k, r, q, v, OptionTypes.EUROPEAN_PUT.value)
            return bsv + a1 * ((s/sstar)**q1)
        else:
            return k - s

    else:

        raise FinError("Phi must equal 1 or -1.")

###############################################################################


@njit(fastmath=True)
def bjerksund_stensland_value(s, t, k, r, q, v, option_type_value):
    """Price American Option using the Bjerksund-Stensland
     approximation (1993) for the Black Scholes Model"""
    if option_type_value == OptionTypes.AMERICAN_CALL.value:
        pass
    elif option_type_value == OptionTypes.AMERICAN_PUT.value:
        # put-call transformation
        s, k, r, q = k, s, r-q, -q
    else:
        return 0.0

    def phi(S, T, gamma, H, X):
        """Eq.(13) in Bjerksund-Stensland approximation (1993)."""
        nonlocal r, q
        lambda0 = (-r + gamma * q + 0.5 * gamma *
                   (gamma - 1.0) * v**2) * T
        d = - (np.log(S/H) + (q + (gamma - 0.5) * v**2) * T) / \
            (v * np.sqrt(t))
        kappa = (2.0 * gamma - 1.0) + (2.0 * q) / v**2
        return (
            np.exp(lambda0) * (S ** gamma)
            * (N(d) - N(d - (2.0 * np.log(X/S)/v/np.sqrt(T))) * ((X/S)**kappa))
        )

    # calc trigger price x_t
    beta = (0.5 - q/(v**2)) + np.sqrt((0.5 - q/(v**2))**2 + 2.0 * r/(v**2))
    # avoid division by zero
    if abs(r-q) < 1.e-10:
        beta = 1.0
        x_t = 1.e10
    else:
        b_infty = k * beta / (beta - 1.0)
        b_0 = max(k, k * r/((r-q)))
        h_t = -(q*t + 2.0 * v * np.sqrt(t)) * (b_0 / (b_infty - b_0))
        x_t = b_0 + (b_infty - b_0) * (1.0 - np.exp(h_t))
    # calc option value
    alpha = (x_t - k) * x_t ** (-beta)
    value = (alpha * (s**beta) - alpha * phi(s, t, beta, x_t, x_t)
             + phi(s, t, 1.0, x_t, x_t) - phi(s, t, 1.0, k, x_t)
             - k * phi(s, t, 0.0, x_t, x_t) + k * phi(s, t, 0.0, k, x_t))

    return value

###############################################################################
