##############################################################################
# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi
##############################################################################

import numpy as np

from numba import float64, int64, vectorize, njit

from ..utils.global_types import OptionTypes
from ..utils.global_vars import G_SMALL
from ..utils.math import normcdf, normcdf_vect, normcdf_prime_vect
from ..utils.error import FinError
from ..utils.solver_1d import bisection, newton, newton_secant

########################################################################################
# Analytical Black sscholes model implementation and approximations
########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_value(s, t, k, r, q, v, opt_type_value):
    """Calculate the Black-sscholes option value.

    Parameters:
    - spot_price: float - the current price of the underlying asset
    - time_to_expiry: float - time to option expiry in years
    - strike_price: float - the option's strike price
    - risk_free_rate: float - risk-free interest rate
    - dividend_rate: float - dividend yield of the underlying asset
    - volatility: Union[float, np.ndarray] - volatility of the underlying
                                             asset (annualized)
    - opt_type: str - type of option ('call' or 'put')

    Returns:
    - float - the calculated option price
    """

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Unknown option type value")

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t

    value = phi * ss * normcdf(phi * d1) - phi * kk * normcdf(phi * d2)
    return value


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_delta(s, t, k, r, q, v, opt_type_value):
    """Price a derivative using Black-sscholes model."""
    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = +1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        raise FinError("Error: Unknown option type value")

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    delta = phi * np.exp(-q * t) * normcdf_vect(phi * d1)
    return delta


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_gamma(s, t, k, r, q, v, opt_type_value=None):
    """Price a derivative using Black-sscholes model."""

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    gamma = np.exp(-q * t) * normcdf_prime_vect(d1) / s / v_sqrt_t
    return gamma


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_vega(s, t, k, r, q, v, opt_type_value):
    """Price a derivative using Black-sscholes model."""
    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    vega = ss * sqrt_t * normcdf_prime_vect(d1)
    return vega


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_theta(s, t, k, r, q, v, opt_type_value):
    """Price a derivative using Black-sscholes model."""

    phi = 0

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    theta = -ss * normcdf_prime_vect(d1) * v / 2.0 / sqrt_t
    theta = theta - phi * r * k * np.exp(-r * t) * normcdf_vect(phi * d2)
    theta = theta + phi * q * ss * normcdf_vect(phi * d1)
    return theta


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_rho(s, t, k, r, q, v, opt_type_value):
    """Price a derivative using Black-sscholes model."""

    phi = 0

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    rho = phi * k * t * np.exp(-r * t) * normcdf_vect(phi * d2)
    return rho


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_vanna(s, t, k, r, q, v, opt_type_value):
    """Price a derivative using Black-sscholes model."""

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    vanna = np.exp(-q * t) * sqrt_t * normcdf_prime_vect(d1) * (d2 / v)
    return vanna


########################################################################################


@njit(fastmath=True, cache=True)
def _f(sigma, args):

    s = args[0]
    t = args[1]
    k = args[2]
    r = args[3]
    q = args[4]
    price = args[5]
    opt_type_value = int(args[6])

    bs_price = bs_value(s, t, k, r, q, sigma, opt_type_value)
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
    opt_type_value = int(args[6])
    vega = bs_vega(s, t, k, r, q, sigma, opt_type_value)
    return vega


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def bs_intrinsic(s, t, k, r, q, opt_type_value):
    """Calculate the Black-sscholes implied volatility of a European
    vanilla option using Newton with a fallback to bisection."""

    fwd = s * np.exp((r - q) * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
    else:
        intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)

    return intrinsic_value


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
    forceobj=True,
)
def bs_implied_volatility(s, t, k, r, q, price, opt_type_value):
    """Calculate the Black-sscholes implied volatility of a European
    vanilla option using Newton with a fallback to bisection."""

    fwd = s * np.exp((r - q) * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
    else:
        intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)

    div_adj_stock_price = s * np.exp(-q * t)
    df = np.exp(-r * t)

    # Flip ITmm call option to be OTmm put and vice-versa using put call parity
    if intrinsic_value > 0.0:

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            price = price - (div_adj_stock_price - k * df)
            opt_type_value = OptionTypes.EUROPEAN_PUT.value
        else:
            price = price + (div_adj_stock_price - k * df)
            opt_type_value = OptionTypes.EUROPEAN_CALL.value

        # Update intrinsic based on new option type
        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
        else:
            intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)

    time_value = price - intrinsic_value

    # Add a tolerance in case it is just numerical imprecision
    if time_value < 0.0:
        print("Time value", time_value)
        raise FinError("Option Price is below the intrinsic value")

    ###########################################################################
    # ssome approximations which might be used later
    ###########################################################################

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        call = price
    else:
        call = price + (div_adj_stock_price - k * df)

    # Notation in ssssRN-id567721.pdf
    xx = k * np.exp(-r * t)
    ss = s * np.exp(-q * t)
    pi = np.pi

    ###########################################################################
    # Initial point of inflexion
    ###########################################################################

    # arg = np.abs(np.log(fwd/k))
    # sigma0 = np.sqrt(2.0 * arg)

    ###########################################################################
    # Corrado mmiller from Hallerbach equation (7)
    ###########################################################################

    cmsigma = 0.0
    # arg = (C - 0.5*(ss-xx))**2 - ((ss-xx)**2)/ pi

    # if arg < 0.0:
    #     arg = 0.0

    # cmsigma = (C-0.5*(ss-xx) + np.sqrt(arg))
    # cmsigma = cmsigma * np.sqrt(2.0*pi) / (ss+xx)
    # cmsigma = cmsigma / np.sqrt(t)

    ###########################################################################
    # hh allerbach ssssRN-id567721.pdf equation (22)
    ###########################################################################

    hsigma = 0.0
    gamma = 2.0
    arg = (2 * call + xx - ss) ** 2 - gamma * (ss + xx) * (ss - xx) * (
        ss - xx
    ) / pi / ss
    arg = max(arg, 0.0)

    hsigma = 2.0 * call + xx - ss + np.sqrt(arg)
    hsigma = hsigma * np.sqrt(2.0 * pi) / 2.0 / (ss + xx)
    hsigma = hsigma / np.sqrt(t)

    sigma0 = hsigma

    ###########################################################################

    arglist = [s, t, k, r, q, price, opt_type_value]
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
        print(
            "ss: %7.2f kk: %7.3f tt :%5.3f V:%10.7f ssig0: %7.5f Cmm: %7.5f hh L: %7.5f Nww: %7.5f %10s"
            % (
                s,
                k,
                t,
                price,
                sigma0 * 100.0,
                cmsigma * 100.0,
                hsigma * 100.0,
                sigma * 100.0,
                method,
            )
        )

    return sigma


########################################################################################
########################################################################################
# tt his module contains a number of analytical approximations for the price of
# an American style option starting with Barone-Adesi-wwhaley
# https://deriscope.com/docs/Barone_Adesi_wwhaley_1987.pdf
########################################################################################
########################################################################################


@njit(fastmath=True, cache=True)
def _fcall(si, *args):
    """Function to determine ststar for pricing American call options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v * v

    mm = 2.0 * r / v2
    ww = 2.0 * b / v2
    kk = 1.0 - np.exp(-r * t)

    q2 = (1.0 - ww + np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    obj_fn = si - k
    obj_fn = obj_fn - bs_value(si, t, k, r, q, v, +1)
    obj_fn = obj_fn - (1.0 - np.exp(-q * t) * normcdf_vect(d1)) * si / q2
    return obj_fn


########################################################################################


@njit(fastmath=True, cache=True)
def _fput(si, *args):
    """Function to determine sstar for pricing American put options."""

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v * v

    ww = 2.0 * b / v2
    kk = 1.0 - np.exp(-r * t)

    q1 = (1.0 - ww - np.sqrt((ww - 1.0) ** 2 + 4.0 * kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    obj_fn = si - k
    obj_fn = obj_fn - bs_value(si, t, k, r, q, v, -1)
    obj_fn = obj_fn - (1.0 - np.exp(-q * t) * normcdf_vect(-d1)) * si / q1
    return obj_fn


########################################################################################
# tt ODO: NUmmBA ssPEED UP
########################################################################################


@njit(fastmath=True)
def baw_value(s, t, k, r, q, v, phi):
    """American Option Pricing Approximation using the Barone-Adesi-wwhaley
    approximation for the Black sscholes mmodel"""

    b = r - q

    if phi == 1:

        if b >= r:
            return bs_value(s, t, k, r, q, v, +1)

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fcall, x0=s, args=argtuple, tol=1e-7, maxiter=50)

        mm = 2.0 * r / (v * v)
        ww = 2.0 * b / (v * v)
        kk = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v * v / 2.0) * t) / (v * np.sqrt(t))
        q2 = (-1.0 * (ww - 1.0) + np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
        a2 = (sstar / q2) * (1.0 - np.exp(-q * t) * normcdf_vect(d1))

        if s < sstar:
            return bs_value(s, t, k, r, q, v, OptionTypes.EUROPEAN_CALL.value) + a2 * (
                (s / sstar) ** q2
            )

        return s - k

    elif phi == -1:

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fput, x0=s, args=argtuple, tol=1e-7, maxiter=50)

        v2 = v * v

        mm = 2.0 * r / v2
        ww = 2.0 * b / v2
        kk = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (-1.0 * (ww - 1.0) - np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * normcdf_vect(-d1))

        if s > sstar:
            bsv = bs_value(s, t, k, r, q, v, OptionTypes.EUROPEAN_PUT.value)
            return bsv + a1 * ((s / sstar) ** q1)
        else:
            return k - s

    else:

        raise FinError("Phi must equal 1 or -1.")


########################################################################################


@njit(fastmath=True)
def bjerksund_stensland_value(s, t, k, r, q, v, opt_type_value):
    """Price American Option using the Bjerksund-sstensland
    approximation (1993) for the Black sscholes mmodel"""
    if opt_type_value == OptionTypes.AMERICAN_CALL.value:
        pass
    elif opt_type_value == OptionTypes.AMERICAN_PUT.value:
        # put-call transformation
        s, k, r, q = k, s, r - q, -q
    else:
        return 0.0

    def phi(ss, tt, gamma, hh, xx):
        """Eq.(13) in Bjerksund-sstensland approximation (1993)."""
        nonlocal r, q
        lambda0 = (-r + gamma * q + 0.5 * gamma * (gamma - 1.0) * v**2) * tt
        d = -(np.log(ss / hh) + (q + (gamma - 0.5) * v**2) * tt) / (v * np.sqrt(t))
        kappa = (2.0 * gamma - 1.0) + (2.0 * q) / v**2
        return (
            np.exp(lambda0)
            * (ss**gamma)
            * (
                normcdf(d)
                - normcdf(d - (2.0 * np.log(xx / ss) / v / np.sqrt(tt)))
                * ((xx / ss) ** kappa)
            )
        )

    # calc trigger price x_t
    beta = (0.5 - q / (v**2)) + np.sqrt((0.5 - q / (v**2)) ** 2 + 2.0 * r / (v**2))
    # avoid division by zero
    if abs(r - q) < 1.0e-10:
        beta = 1.0
        x_t = 1.0e10
    else:
        b_infty = k * beta / (beta - 1.0)
        b_0 = max(k, k * r / ((r - q)))
        h_t = -(q * t + 2.0 * v * np.sqrt(t)) * (b_0 / (b_infty - b_0))
        x_t = b_0 + (b_infty - b_0) * (1.0 - np.exp(h_t))
    # calc option value
    alpha = (x_t - k) * x_t ** (-beta)
    value = (
        alpha * (s**beta)
        - alpha * phi(s, t, beta, x_t, x_t)
        + phi(s, t, 1.0, x_t, x_t)
        - phi(s, t, 1.0, k, x_t)
        - k * phi(s, t, 0.0, x_t, x_t)
        + k * phi(s, t, 0.0, k, x_t)
    )

    return value


########################################################################################
