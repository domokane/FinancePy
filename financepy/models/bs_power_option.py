# Copyright (C) 2020 Dominic O'Kane, G Poorna Prudhvi

import numpy as np
from numba import njit

from ..utils.global_types import OptionTypes
from ..utils.global_vars import G_SMALL
from ..utils.math import normcdf

########################################################################################


@njit(fastmath=True, cache=True)
def power_option_value(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
    n: float,
) -> float:
    """
    Price a European power option with payoff:

        Call: max(S_T**n - K, 0)
        Put:  max(K - S_T**n, 0)

    under Black-Scholes.
    """

    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        return np.nan

    if n <= 0.0:
        return np.nan

    s = max(s, G_SMALL)
    k = max(k, G_SMALL)
    t = max(t, G_SMALL)
    v = max(v, G_SMALL)

    sqrt_t = np.sqrt(t)

    # Critical stock level:
    # S_T**n > K  <=>  S_T > K**(1/n)
    k_root = k ** (1.0 / n)

    # Distribution threshold
    d2 = (
        np.log(s / k_root)
        + (r - q - 0.5 * v * v) * t
    ) / (v * sqrt_t)

    # Shift caused by the S_T**n moment
    d1 = d2 + n * v * sqrt_t

    # Discounted n-th moment of S_T
    stock_term = (
        s ** n
        * np.exp(
            (
                (n - 1.0) * r
                - n * q
                + 0.5 * n * (n - 1.0) * v * v
            ) * t
        )
    )

    strike_term = k * np.exp(-r * t)

    return (
        phi * stock_term * normcdf(phi * d1)
        - phi * strike_term * normcdf(phi * d2)
    )


@njit(fastmath=True, cache=True)
def power_option_delta(
    s: float,
    t: float,
    k: float,
    r: float,
    q: float,
    v: float,
    opt_type_value: int,
    n: float,
) -> float:
    """Analytical spot delta of the existing European power-option price.

    For positive power n, delta = phi * n * M_n / S * N(phi * d1),
    where M_n = S**n * exp(((n-1)*r - n*q + n*(n-1)*v*v/2)*t)
    is the discounted n-th stock-price moment. Uses the same small-
    time/volatility conventions as power_option_value.
    """
    if opt_type_value == OptionTypes.EUROPEAN_CALL.value:
        phi = 1.0
    elif opt_type_value == OptionTypes.EUROPEAN_PUT.value:
        phi = -1.0
    else:
        return np.nan

    if n <= 0.0:
        return np.nan

    s = max(s, G_SMALL)
    k = max(k, G_SMALL)
    t = max(t, G_SMALL)
    v = max(v, G_SMALL)

    sqrt_t = np.sqrt(t)
    k_root = k ** (1.0 / n)
    d1 = (
        (np.log(s / k_root) + (r - q - 0.5 * v * v) * t)
        / (v * sqrt_t)
        + n * v * sqrt_t
    )

    stock_term = s ** n * np.exp(
        ((n - 1.0) * r - n * q
         + 0.5 * n * (n - 1.0) * v * v) * t
    )

    return phi * (n * stock_term / s) * normcdf(phi * d1)

#########################################################################################
