# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from enum import Enum

import numpy as np
from numba import njit, float64

from ..utils.math import normcdf
from ..utils.error import FinError

# Parametric functions for option volatility to use in a Black-Scholes model

########################################################################################


class VolFuncTypes(Enum):

    CLARK = 0
    SABR = 1
    SABR_BETA_ONE = 2
    SABR_BETA_HALF = 3
    BBG = 4
    CLARK5 = 5
    SVI = 6
    SSVI = 7


########################################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def vol_function_clark(
    params: np.ndarray,
    f: float,
    k: float,
    t: float
) -> float:
    """Volatility Function in book by Iain Clark generalised to allow for
    higher than quadratic power. Care needs to be taken to avoid overfitting.
    The exact reference is Clark Page 59."""

    if f < 0.0:
        print("f:", f)
        raise FinError("Forward is negative")

    if k < 0.0:
        print("k:", k)
        raise FinError("Strike is negative")

    x = np.log(f / k)
    sigma0 = np.exp(params[0])
    arg = x / (sigma0 * np.sqrt(t))
    deltax = normcdf(arg) - 0.50  # The -0.50 seems to be missing in book

    f = 0.0
    for i in range(0, len(params)):
        f += params[i] * (deltax**i)

    return np.exp(f)


########################################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def vol_function_bloomberg(
    params: np.ndarray,
    f: float,
    k: float,
    t: float
) -> float:
    """Volatility Function similar to the one used by Bloomberg. It is
    a quadratic function in the spot delta of the option. It can therefore
    go negative so it requires a good initial guess when performing the
    fitting to avoid this happening. The first parameter is the quadratic
    coefficient i.e. sigma(K) = a * D * D + b * D + c where a = params[0],
    b = params[1], c = params[2] and D is the spot delta."""

    num_params = len(params)

    # Rather than pass in the ATM vol, I imply it from the delta=0.50 curve
    sigma = 0.0
    for i in range(0, len(params)):
        pwr = num_params - i - 1
        sigma += params[i] * ((0.50) ** pwr)

    vsqrtt = sigma * np.sqrt(t)

    d1 = np.log(f / k) / vsqrtt + vsqrtt / 2.0
    delta = normcdf(d1)

    v = 0.0
    for i in range(0, len(params)):
        pwr = num_params - i - 1
        v += params[i] * (delta**pwr)

    return v


# I do not jit this so it can be called from a notebook with a vector of strike
# Also, if I vectorise it it fails as it cannot handle a numpy array as input

########################################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def vol_function_svi(
    params: np.ndarray,
    f: float,
    k: float,
    t: float
) -> float:
    """Volatility Function proposed by Gatheral in 2004. Increasing a results
    in a vertical translation of the smile in the positive direction.
    Increasing b decreases the angle between the put and call wing, i.e.
    tightens the smile. Increasing rho results in a counter-clockwise rotation
    of the smile. Increasing m results in a horizontal translation of the smile
    in the positive direction. Increasing sigma reduces the at-the-money
    curvature of the smile."""

    x = np.log(f / k)

    a = params[0]
    b = params[1]
    rho = params[2]
    m = params[3]
    sigma = params[4]

    vart = a + b * (rho * (x - m) + np.sqrt((x - m) ** 2 + sigma * sigma))
    v = np.sqrt(vart / t)
    return v


# Gatheral ssvi surface svi and equivalent local volatility
# Code from https://wwwf.imperial.ac.uk/~ajacquie/IC_AMDP/IC_AMDP_Docs/Code/ssvi.pdf

########################################################################################


@njit(float64(float64, float64), fastmath=True, cache=True)
def phi_ssvi(
    theta: float,
    gamma: float
) -> float:

    if abs(gamma) < 1e-8:
        gamma = 1e-8

    if abs(theta) < 1e-8:
        theta = 1e-8

    phi = (1.0 / gamma / theta) * (1.0 - (1.0 - np.exp(-gamma * theta)) / gamma / theta)
    return phi


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def ssvi(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:
    """This is the total variance w = sigma(t) x sigma(t) (0,t) x t"""

    theta = sigma * sigma * t
    p = phi_ssvi(theta, gamma)
    px = p * x
    gg = px + rho
    v = 0.5 * theta * (1.0 + rho * px + np.sqrt(gg**2 + 1.0 - rho * rho))
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def ssvi1(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    # First derivative with respect to x
    theta = sigma * sigma * t
    p = phi_ssvi(theta, gamma)
    px = p * x
    v = 0.5 * theta * p * (px + rho * np.sqrt(px**2 + 2.0 * px * rho + 1.0) + rho)
    v = v / np.sqrt(px**2 + 2.0 * px * rho + 1.0)
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def ssvi2(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    # Second derivative with respect to x
    theta = sigma * sigma * t
    p = phi_ssvi(theta, gamma)
    px = p * x
    v = 0.5 * theta * p * p * (1.0 - rho * rho)
    v = v / ((px**2 + 2.0 * px * rho + 1.0) * np.sqrt(px**2 + 2.0 * px * rho + 1.0))
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def ssvit(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    # First derivative with respect to t, by central difference
    eps = 0.0001
    ssvitplus = ssvi(x, gamma, sigma, rho, t + eps)
    ssvitminus = ssvi(x, gamma, sigma, rho, t - eps)
    deriv = (ssvitplus - ssvitminus) / 2.0 / eps
    return deriv


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def g(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    w = ssvi(x, gamma, sigma, rho, t)

    if abs(w) < 1e-10:
        w = 1e-10

    w1 = ssvi1(x, gamma, sigma, rho, t)
    w2 = ssvi2(x, gamma, sigma, rho, t)
    xwv = x * w1 / w
    v = (1.0 - 0.5 * xwv) ** 2 - 0.25 * w1 * w1 * (0.25 + 1.0 / w) + 0.5 * w2
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def dminus(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    vsqrt = np.sqrt(ssvi(x, gamma, sigma, rho, t))
    v = -x / vsqrt - 0.5 * vsqrt
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def density_ssvi(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    dm = dminus(x, gamma, sigma, rho, t)
    v = g(x, gamma, sigma, rho, t) * np.exp(-0.5 * dm * dm)
    v = v / np.sqrt(2.0 * np.pi * ssvi(x, gamma, sigma, rho, t))
    return v


########################################################################################


@njit(
    float64(float64, float64, float64, float64, float64),
    fastmath=True,
    cache=True,
)
def ssvi_local_varg(
    x: float,
    gamma: float,
    sigma: float,
    rho: float,
    t: float
) -> float:

    # Compute the equivalent ssvi local variance
    num = ssvit(x, gamma, sigma, rho, t)
    den = g(x, gamma, sigma, rho, t)
    var = num / den
    return var


########################################################################################


@njit(float64(float64[:], float64, float64, float64), fastmath=True, cache=True)
def vol_function_ssvi(
    params: np.ndarray,
    f: float,
    k: float,
    t: float
) -> float:
    """Volatility Function proposed by Gatheral in 2004."""

    gamma = params[0]
    sigma = params[1]
    rho = params[2]

    x = np.log(f / k)

    vart = ssvi_local_varg(x, gamma, sigma, rho, t)
    vart = max(vart, 0.0)

    sigma = np.sqrt(vart)

    return sigma
