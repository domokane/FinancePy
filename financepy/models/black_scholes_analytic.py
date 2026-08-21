# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi

import numpy as np

from numba import float64, int64, vectorize, njit

from ..utils.global_types import OptionTypes
from ..utils.global_vars import G_SMALL
from ..utils.math import normcdf, normcdf_vect, normcdf_prime_vect
from ..utils.solver_1d import bisection, newton, newton_secant

# Analytical Black-Scholes model implementation and approximations for American style

########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64)],
    fastmath=True,
    cache=True,
)
def fwd(
    s: float,
    t: float,
    r: float,
    q: float,
) -> float:
    """Return the forward price.
    Parameters:
    - spot_price: float - the current price of the underlying asset
    - time_to_expiry: float - time to option expiry in years
    - risk_free_rate: float - risk-free interest rate
    - dividend_rate: float - dividend yield of the underlying asset

    Returns:
    - float - the calculated fwd price
    """

    fwd = s * np.exp((r - q) * t)
    return fwd


#########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def european_value(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price a european style derivative using Black-Scholes model.
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
        return np.nan

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t

    v = phi * ss * normcdf(phi * d1) - phi * kk * normcdf(phi * d2)
    return v


#########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)],
    fastmath=True,
    cache=True,
)
def d1(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
) -> float:
    """Return d1 and d2
    - spot_price: float - the current price of the underlying asset
    - time_to_expiry: float - time to option expiry in years
    - strike_price: float - the option's strike price
    - risk_free_rate: float - risk-free interest rate
    - dividend_rate: float - dividend yield of the underlying asset
    - volatility: Union[float, np.ndarray] - volatility of the underlying
                                             asset (annualized)

    Returns:
    - float - the calculated d1
    """

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0

    return d1


#########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64)],
    fastmath=True,
    cache=True,
)
def d2(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
) -> float:
    """Return d1 and d2
    - spot_price: float - the current price of the underlying asset
    - time_to_expiry: float - time to option expiry in years
    - strike_price: float - the option's strike price
    - risk_free_rate: float - risk-free interest rate
    - dividend_rate: float - dividend yield of the underlying asset
    - volatility: Union[float, np.ndarray] - volatility of the underlying
                                             asset (annualized)
    Returns:
    - float - the calculated d2
    """

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t

    return d2


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def delta(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price a derivative using Black-Scholes model."""
    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = +1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        return np.nan

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d = phi * np.exp(-q * t) * normcdf_vect(phi * d1)
    return d


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def gamma(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int = None,
) -> float:
    """Price a derivative using Black-sscholes model."""

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    v_sqrt_t = v * np.sqrt(t)
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    g = np.exp(-q * t) * normcdf_prime_vect(d1) / s / v_sqrt_t
    return g


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def vega(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price a derivative using Black-Scholes model."""
    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    v = ss * sqrt_t * normcdf_prime_vect(d1)
    return v


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def theta(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price a derivative using Black-Scholes model."""

    phi = 0

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        return np.nan

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

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
def rho(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price a derivative using Black-Scholes model."""

    phi = 0

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        return np.nan

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    r = phi * k * t * np.exp(-r * t) * normcdf_vect(phi * d2)
    return r


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def vanna(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Calculate European option vanna under the Black–Scholes model."""

    k = np.maximum(k, G_SMALL)
    t = np.maximum(t, G_SMALL)
    v = np.maximum(v, G_SMALL)
    s = np.maximum(s, G_SMALL)

    sqrt_t = np.sqrt(t)
    v_sqrt_t = v * sqrt_t
    ss = s * np.exp(-q * t)
    kk = k * np.exp(-r * t)
    d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
    d2 = d1 - v_sqrt_t
    vanna = -np.exp(-q * t) * normcdf_prime_vect(d1) * (d2 / v)
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

    bs_price = european_value(s, t, k, r, q, sigma, opt_type_value)
    obj = bs_price - price
    return obj


########################################################################################


@njit(fastmath=True, cache=True)
def _fvega(sigma, args):

    s = args[0]
    t = args[1]
    k = args[2]
    r = args[3]
    q = args[4]
    opt_type_value = int(args[6])
    v = vega(s, t, k, r, q, sigma, opt_type_value)
    return v


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def intrinsic(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    opt_type_value: int,
) -> float:
    """Calculate the Black-Scholes implied volatility of a European
    vanilla option using Newton with a fallback to bisection."""

    fwd = s * np.exp((r - q) * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)
    else:
        intrinsic_value = np.nan

    return intrinsic_value


########################################################################################


# @vectorize(
#     [float64(float64, float64, float64, float64, float64, float64, int64)],
#     fastmath=True,
#     cache=True,
# )
def implied_volatility(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    price: float,
    opt_type_value: int,
) -> float:
    """Calculate the Black-Scholes implied volatility of a European
    vanilla option using Newton with a fallback to bisection."""

    if t <= 0.0:
        return np.nan

    fwd = s * np.exp((r - q) * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)
    else:
        return np.nan

    div_adj_stock_price = s * np.exp(-q * t)
    df = np.exp(-r * t)

    # Flip ITmm call option to be OTmm put and vice-versa using put call parity
    if intrinsic_value > 0.0:

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            price = price - (div_adj_stock_price - k * df)
            opt_type_value = OptionTypes.EUROPEAN_PUT.value
        elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
            price = price + (div_adj_stock_price - k * df)
            opt_type_value = OptionTypes.EUROPEAN_CALL.value
        else:
            return np.nan

        # Update intrinsic based on new option type
        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
        elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
            intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)
        else:
            return np.nan

    time_value = price - intrinsic_value

    # Add a tolerance in case it is just numerical imprecision
    if time_value < -1.0 * G_SMALL:
        return np.nan

    time_value = max(time_value, 0.0)

    # some approximations which might be used later

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        call = price
    else:
        call = price + (div_adj_stock_price - k * df)

    # Notation in SSRN-id567721.pdf
    xx = k * np.exp(-r * t)
    ss = s * np.exp(-q * t)
    pi = np.pi

    # Initial point of inflexion

    # arg = np.abs(np.log(fwd/k))
    # sigma0 = np.sqrt(2.0 * arg)

    # Corrado mmiller from Hallerbach equation (7)

    cmsigma = 0.0
    # arg = (C - 0.5*(ss-xx))**2 - ((ss-xx)**2)/ pi

    # if arg < 0.0:
    #     arg = 0.0

    # cmsigma = (C-0.5*(ss-xx) + np.sqrt(arg))
    # cmsigma = cmsigma * np.sqrt(2.0*pi) / (ss+xx)
    # cmsigma = cmsigma / np.sqrt(t)

    # hh allerbach ssssRN-id567721.pdf equation (22)

    hsigma = 0.0
    gam = 2.0
    arg = (2 * call + xx - ss) ** 2 - gam * (ss + xx) * (ss - xx) * (ss - xx) / pi / ss
    arg = max(arg, 0.0)

    hsigma = 2.0 * call + xx - ss + np.sqrt(arg)
    hsigma = hsigma * np.sqrt(2.0 * pi) / 2.0 / (ss + xx)
    hsigma = hsigma / np.sqrt(t)

    sigma0 = hsigma

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

    debug = False
    if debug:
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


# This module contains a number of analytical approximations for the price of
# an American style option starting with Barone-Adesi-Whaley
# https://deriscope.com/docs/Barone_Adesi_Whaley_1987.pdf

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
    obj_fn = obj_fn - european_value(si, t, k, r, q, v, OptionTypes.EUROPEAN_CALL.value)
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

    mm = 2.0 * r / v2
    ww = 2.0 * b / v2
    kk = 1.0 - np.exp(-r * t)

    q1 = (1.0 - ww - np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
    obj_fn = si - k
    obj_fn = obj_fn + european_value(
        si,
        t,
        k,
        r,
        q,
        v,
        OptionTypes.EUROPEAN_PUT.value,
    )
    obj_fn = obj_fn - (1.0 - np.exp(-q * t) * normcdf_vect(-d1)) * si / q1
    return obj_fn


########################################################################################


@njit(fastmath=True)
def baw_value(s, t, k, r, q, v, opt_type_value):
    """American Option Pricing Approximation using the Barone-Adesi-wwhaley
    approximation for the Black-Scholes mmodel"""

    b = r - q

    if opt_type_value == OptionTypes.AMERICAN_CALL.value:

        if t <= G_SMALL:
            return max(s - k, 0.0)

        euro_type = OptionTypes.EUROPEAN_CALL.value

        if b >= r:
            return european_value(s, t, k, r, q, v, euro_type)

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fcall, x0=s, args=argtuple, tol=1e-7, maxiter=50)

        mm = 2.0 * r / (v * v)
        ww = 2.0 * b / (v * v)
        kk = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v * v / 2.0) * t) / (v * np.sqrt(t))
        q2 = (-1.0 * (ww - 1.0) + np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
        a2 = (sstar / q2) * (1.0 - np.exp(-q * t) * normcdf_vect(d1))

        if s < sstar:
            return european_value(s, t, k, r, q, v, euro_type) + a2 * ((s / sstar) ** q2)

        return s - k

    elif opt_type_value == OptionTypes.AMERICAN_PUT.value:

        if t <= G_SMALL:
            return max(k - s, 0.0)

        euro_type = OptionTypes.EUROPEAN_PUT.value

        argtuple = (t, k, r, q, v)

        sstar = newton_secant(_fput, x0=k, args=argtuple, tol=1e-7, maxiter=50)

        v2 = v * v

        mm = 2.0 * r / v2
        ww = 2.0 * b / v2
        kk = 1.0 - np.exp(-r * t)
        d1 = (np.log(sstar / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))
        q1 = (-1.0 * (ww - 1.0) - np.sqrt((ww - 1.0) ** 2 + 4.0 * mm / kk)) / 2.0
        a1 = -(sstar / q1) * (1 - np.exp(-q * t) * normcdf_vect(-d1))

        if s > sstar:
            bsv = european_value(s, t, k, r, q, v, euro_type)
            return bsv + a1 * ((s / sstar) ** q1)
        else:
            return k - s

    return np.nan


########################################################################################


@njit(fastmath=True)
def bjerksund_stensland_value(s, t, k, r, q, v, opt_type_value):
    """Price American Option using the Bjerksund-Stensland
    approximation (1993) for the Black-Scholes mmodel"""
    if opt_type_value == OptionTypes.AMERICAN_CALL.value:
        pass
    elif opt_type_value == OptionTypes.AMERICAN_PUT.value:
        # put-call transformation
        s, k, r, q = k, s, q, r
    else:
        return np.nan

    b = r - q

    ####################################################################################

    def phi(ss, tt, gamma, hh, xx):
        """Eq.(13) in Bjerksund-Stensland approximation (1993)."""
        nonlocal r, b
        v2 = v * v
        sqrt_t = np.sqrt(t)

        lambda0 = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * v2) * tt
        d = -(np.log(ss / hh) + (b + (gamma - 0.5) * v2) * tt) / (v * sqrt_t)
        kappa = (2.0 * gamma - 1.0) + (2.0 * b) / v2
        return (
            np.exp(lambda0)
            * (ss**gamma)
            * (normcdf(d) - normcdf(d - (2.0 * np.log(xx / ss) / v / sqrt_t)) * ((xx / ss) ** kappa))
        )

    # calc trigger price x_t
    beta = (0.5 - b / (v**2)) + np.sqrt((b / (v**2) - 0.5) ** 2 + 2.0 * r / v**2)
    # avoid division by zero
    if abs(b) < 1.0e-10:
        beta = 1.0
        x_t = 1.0e10
    else:
        b_infty = k * beta / (beta - 1.0)
        b_0 = max(k, k * r / b)
        h_t = -(b * t + 2.0 * v * np.sqrt(t)) * (b_0 / (b_infty - b_0))
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


###############################################################################


@njit(cache=True)
def _bs93_phi(s, t, gamma, h, i, r, b, v):
    """Equation (13) in the Bjerksund–Stensland 1993 approximation."""

    v2 = v * v
    sqrt_t = np.sqrt(t)

    lambda_ = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * v2) * t

    d = -(np.log(s / h) + (b + (gamma - 0.5) * v2) * t) / (v * sqrt_t)

    kappa = 2.0 * gamma - 1.0 + 2.0 * b / v2

    reflected_d = d - 2.0 * np.log(i / s) / (v * sqrt_t)

    return np.exp(lambda_) * s**gamma * (normcdf(d) - (i / s) ** kappa * normcdf(reflected_d))


########################################################################################


@njit(cache=True)
def bjerksund_stensland_value2(
    s,
    t,
    k,
    r,
    q,
    v,
    opt_type_value,
):
    """Bjerksund–Stensland 1993 American-option approximation.

    This implementation assumes:
        s > 0
        k > 0
        t >= 0
        v > 0
        r >= 0
        q >= 0
    """

    is_call = opt_type_value == OptionTypes.AMERICAN_CALL.value
    is_put = opt_type_value == OptionTypes.AMERICAN_PUT.value

    if not is_call and not is_put:
        return np.nan

    # Handle expiry before performing the put-call transformation.
    if t <= G_SMALL:
        if is_call:
            return max(s - k, 0.0)

        return max(k - s, 0.0)

    if s <= 0.0 or k <= 0.0 or v <= 0.0:
        return np.nan

    # This simple implementation does not handle the negative-rate
    # double-boundary regimes.
    if r < 0.0 or q < 0.0:
        return np.nan

    if is_put:
        # American put-call transformation:
        #
        # P(S, K, r, q) = C(K, S, q, r)
        transformed_s = k
        transformed_k = s
        transformed_r = q
        transformed_q = r

        s = transformed_s
        k = transformed_k
        r = transformed_r
        q = transformed_q

    b = r - q
    v2 = v * v
    sqrt_t = np.sqrt(t)

    euro_call_type = OptionTypes.EUROPEAN_CALL.value

    # After the put-call transformation we are always pricing a call.
    #
    # For q <= 0, early exercise of the call is not optimal in the
    # standard non-negative-rate setting.
    if q <= 0.0:
        return european_value(
            s,
            t,
            k,
            r,
            q,
            v,
            euro_call_type,
        )

    beta = 0.5 - b / v2 + np.sqrt((b / v2 - 0.5) ** 2 + 2.0 * r / v2)

    if not np.isfinite(beta) or beta <= 1.0:
        return np.nan

    b_infinity = k * beta / (beta - 1.0)

    # Correct formula:
    #
    # B0 = max(K, rK / (r-b))
    #
    # and r-b = q.
    b_0 = max(k, k * r / q)

    denominator = b_infinity - b_0

    if denominator <= 0.0 or not np.isfinite(denominator):
        return np.nan

    h_t = -(b * t + 2.0 * v * sqrt_t) * b_0 / denominator

    x_t = b_0 + denominator * (-np.expm1(h_t))

    if not np.isfinite(x_t) or x_t <= k:
        return np.nan

    # The continuation formula is only valid below the trigger.
    if s >= x_t:
        return s - k

    alpha = (x_t - k) * x_t ** (-beta)

    value = (
        alpha * s**beta
        - alpha * _bs93_phi(s, t, beta, x_t, x_t, r, b, v)
        + _bs93_phi(s, t, 1.0, x_t, x_t, r, b, v)
        - _bs93_phi(s, t, 1.0, k, x_t, r, b, v)
        - k * _bs93_phi(s, t, 0.0, x_t, x_t, r, b, v)
        + k * _bs93_phi(s, t, 0.0, k, x_t, r, b, v)
    )

    european_value_transformed = european_value(
        s,
        t,
        k,
        r,
        q,
        v,
        euro_call_type,
    )

    intrinsic_value = max(s - k, 0.0)

    # Numerical safeguards. Under the transformation, these are also
    # the correct European and intrinsic lower bounds for the put.
    return max(
        value,
        european_value_transformed,
        intrinsic_value,
    )


###############################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=False,
)
def value(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
) -> float:
    """Price European or American vanilla options."""

    BAW = False  # True
    #    BAW = True

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        return european_value(s, t, k, r, q, v, opt_type_value)

    if opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        return european_value(s, t, k, r, q, v, opt_type_value)

    if opt_type_value == OptionTypes.AMERICAN_CALL.value:
        if BAW:
            return baw_value(s, t, k, r, q, v, opt_type_value)
        else:
            return bjerksund_stensland_value2(s, t, k, r, q, v, opt_type_value)

    if opt_type_value == OptionTypes.AMERICAN_PUT.value:
        if BAW:
            return baw_value(s, t, k, r, q, v, opt_type_value)
        else:
            return bjerksund_stensland_value2(s, t, k, r, q, v, opt_type_value)

    return np.nan


########################################################################################
