# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi

import numpy as np

from numba import float64, int64, vectorize, njit

from ..utils.global_types import OptionTypes
from ..utils.global_vars import G_SMALL
from ..utils.math import normcdf, normcdf_vect, normcdf_prime_vect

# Analytical Black-Scholes model implementation and approximations for American style

########################################################################################


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


# @vectorize(
#     [float64(float64, float64, float64, float64, float64, float64)],
#     fastmath=True,
#     cache=True,
# )
# def d1(
#     s: float,
#     t: float,
#     k: float,
#     r: float,
#     q: float,
#     v: float,
# ) -> float:
#     """Return d1 and d2
#     - spot_price: float - the current price of the underlying asset
#     - time_to_expiry: float - time to option expiry in years
#     - strike_price: float - the option's strike price
#     - risk_free_rate: float - risk-free interest rate
#     - dividend_rate: float - dividend yield of the underlying asset
#     - volatility: Union[float, np.ndarray] - volatility of the underlying
#                                              asset (annualized)

#     Returns:
#     - float - the calculated d1
#     """

#     k = np.maximum(k, G_SMALL)
#     t = np.maximum(t, G_SMALL)
#     v = np.maximum(v, G_SMALL)
#     s = np.maximum(s, G_SMALL)

#     v_sqrt_t = v * np.sqrt(t)
#     ss = s * np.exp(-q * t)
#     kk = k * np.exp(-r * t)
#     d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
#     return d1


# #########################################################################################


# @vectorize(
#     [float64(float64, float64, float64, float64, float64, float64)],
#     fastmath=True,
#     cache=True,
# )
# def d2(
#     s: float,
#     t: float,
#     k: float,
#     r: float,
#     q: float,
#     v: float,
# ) -> float:
#     """Return d1 and d2
#     - spot_price: float - the current price of the underlying asset
#     - time_to_expiry: float - time to option expiry in years
#     - strike_price: float - the option's strike price
#     - risk_free_rate: float - risk-free interest rate
#     - dividend_rate: float - dividend yield of the underlying asset
#     - volatility: Union[float, np.ndarray] - volatility of the underlying
#                                              asset (annualized)
#     Returns:
#     - float - the calculated d2
#     """

#     k = np.maximum(k, G_SMALL)
#     t = np.maximum(t, G_SMALL)
#     v = np.maximum(v, G_SMALL)
#     s = np.maximum(s, G_SMALL)

#     v_sqrt_t = v * np.sqrt(t)
#     ss = s * np.exp(-q * t)
#     kk = k * np.exp(-r * t)
#     d1 = np.log(ss / kk) / v_sqrt_t + v_sqrt_t / 2.0
#     d2 = d1 - v_sqrt_t

#     return d2


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


# @njit(fastmath=True, cache=True)
# def _f(sigma, args):

#     s = args[0]
#     t = args[1]
#     k = args[2]
#     r = args[3]
#     q = args[4]
#     price = args[5]
#     opt_type_value = int(args[6])

#     bs_price = european_value(s, t, k, r, q, sigma, opt_type_value)
#     obj = bs_price - price
#     return obj

########################################################################################


@njit(fastmath=True, cache=True)
def _f_iv(sigma, s, t, k, r, q, price, opt_type_value):
    """Black-Scholes implied-volatility objective function."""

    bs_price = european_value(
        s,
        t,
        k,
        r,
        q,
        sigma,
        opt_type_value,
    )

    return bs_price - price


########################################################################################


@njit(fastmath=True, cache=True)
def _secant_iv(
    x0,
    s,
    t,
    k,
    r,
    q,
    price,
    opt_type_value,
    vol_tol=1.0e-8,
    price_tol=1.0e-10,
    maxiter=50,
):
    """Secant solver specialized for Black-Scholes implied volatility."""

    p0 = x0

    eps = 1.0e-4
    p1 = x0 * (1.0 + eps)
    p1 += eps if p1 >= 0.0 else -eps

    q0 = _f_iv(
        p0, s, t, k, r, q, price, opt_type_value
    )
    q1 = _f_iv(
        p1, s, t, k, r, q, price, opt_type_value
    )

    if abs(q0) <= price_tol:
        return p0

    if abs(q1) <= price_tol:
        return p1

    if abs(q1) < abs(q0):
        p0, p1 = p1, p0
        q0, q1 = q1, q0

    for _ in range(maxiter):

        if q1 == q0:
            return np.nan

        if abs(q1) > abs(q0):
            ratio = q0 / q1
            p = (-ratio * p1 + p0) / (1.0 - ratio)
        else:
            ratio = q1 / q0
            p = (-ratio * p0 + p1) / (1.0 - ratio)

        if not np.isfinite(p) or p <= 0.0:
            return np.nan

        qp = _f_iv(
            p, s, t, k, r, q, price, opt_type_value
        )

        # Check the pricing equation as well as the volatility step.
        if abs(qp) <= price_tol:
            return p

        if abs(p - p1) <= vol_tol:
            return p

        p0 = p1
        q0 = q1
        p1 = p
        q1 = qp

    return np.nan


########################################################################################


@njit(fastmath=True, cache=True)
def _bisection_iv(
    x1,
    x2,
    s,
    t,
    k,
    r,
    q,
    price,
    opt_type_value,
    vol_tol=1.0e-8,
    price_tol=1.0e-10,
    maxiter=100,
):
    """Bisection fallback for Black-Scholes implied volatility."""

    f1 = _f_iv(
        x1, s, t, k, r, q, price, opt_type_value
    )
    f2 = _f_iv(
        x2, s, t, k, r, q, price, opt_type_value
    )

    if abs(f1) <= price_tol:
        return x1

    if abs(f2) <= price_tol:
        return x2

    # No root is bracketed.
    if f1 * f2 > 0.0:
        return np.nan

    for _ in range(maxiter):

        xmid = 0.5 * (x1 + x2)

        fmid = _f_iv(
            xmid, s, t, k, r, q, price, opt_type_value
        )

        if abs(fmid) <= price_tol:
            return xmid

        if abs(x2 - x1) <= vol_tol:
            return xmid

        if f1 * fmid < 0.0:
            x2 = xmid
            f2 = fmid
        else:
            x1 = xmid
            f1 = fmid

    return 0.5 * (x1 + x2)


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
    """Return the discounted intrinsic value of a European option."""

    fwd = s * np.exp((r - q) * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        intrinsic_value = np.exp(-r * t) * max(fwd - k, 0.0)
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        intrinsic_value = np.exp(-r * t) * max(k - fwd, 0.0)
    else:
        intrinsic_value = np.nan

    return intrinsic_value


########################################################################################


@vectorize(
    [float64(float64, float64, float64, float64, float64, float64, int64)],
    fastmath=True,
    cache=True,
)
def implied_volatility(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    price: float,
    opt_type_value: int,
) -> float:
    """Calculate Black-Scholes implied volatility for a European option.

    Uses a secant root solver initialized with the Hallerbach
    approximation and falls back to bisection if required.
    """

    if t <= 0.0 or s <= 0.0 or k <= 0.0:
        return np.nan

    xx = k * np.exp(-r * t)
    ss = s * np.exp(-q * t)

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        lower_bound = max(ss - xx, 0.0)
        upper_bound = ss

    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        lower_bound = max(xx - ss, 0.0)
        upper_bound = xx

    else:
        return np.nan

    # Allow a small amount of floating-point noise around the
    # no-arbitrage price bounds.
    vol_tol = 1.0e-6
    price_tol = 1.0e-10

    if price < lower_bound - price_tol:
        return np.nan

    if price > upper_bound + price_tol:
        return np.nan

    # Clamp tiny numerical violations back onto the valid interval.
    price = min(max(price, lower_bound), upper_bound)

    if abs(price - lower_bound) <= price_tol:
        return 0.0

    if abs(price - upper_bound) <= price_tol:
        return np.nan

    intrinsic_value = lower_bound

    # Convert an ITM option into its OTM parity equivalent.
    # Implied volatility is unchanged by put-call parity.
    if intrinsic_value > 0.0:

        parity = ss - xx

        if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
            price = price - parity
            opt_type_value = OptionTypes.EUROPEAN_PUT.value
        else:
            price = price + parity
            opt_type_value = OptionTypes.EUROPEAN_CALL.value

    # Convert the option price to its equivalent call price.
    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        call = price
    else:
        call = price + ss - xx

    # Hallerbach approximation for the initial volatility estimate.
    gam = 2.0

    arg = (
        (2.0 * call + xx - ss) ** 2
        - gam * (ss + xx) * (ss - xx) ** 2
        / (np.pi * ss)
    )

    arg = max(arg, 0.0)

    sigma0 = (
        (2.0 * call + xx - ss + np.sqrt(arg))
        * np.sqrt(2.0 * np.pi)
        / (2.0 * (ss + xx) * np.sqrt(t))
    )

    if not np.isfinite(sigma0) or sigma0 <= 0.0:
        sigma0 = 0.20

    sigma = _secant_iv(
        sigma0,
        s,
        t,
        k,
        r,
        q,
        price,
        opt_type_value,
        vol_tol=vol_tol,
        price_tol=price_tol,
        maxiter=50,
    )

    if not np.isfinite(sigma) or sigma <= 0.0:
        sigma = _bisection_iv(
            1.0e-4,
            10.0,
            s,
            t,
            k,
            r,
            q,
            price,
            opt_type_value,
            vol_tol=vol_tol,
            price_tol=price_tol,
            maxiter=100,
        )

    return sigma


########################################################################################
