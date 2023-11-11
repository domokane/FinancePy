##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import float64, int64, vectorize

from ..utils.error import FinError
from ..utils.global_types import EquityBarrierTypes
from ..utils.math import N


# Calculates the Barrier option price using an Analytical Approach
# and the Black Scholes Model

@vectorize([float64(float64, float64, float64, float64, float64, float64,
                    float64, int64, int64)], fastmath=True, cache=True)
def value_bs(t, k, h, s, r, q, v, option_type, nobs):
    """ This values a single option. Because of its structure it cannot
    easily be vectorised which is why it has been wrapped.
    # number of observations per year
    """
    lnS0k = np.log(s / k)
    sqrtT = np.sqrt(t)

    sigma_rt_t = v * sqrtT
    v2 = v * v
    mu = r - q
    d1 = (lnS0k + (mu + v2 / 2.0) * t) / sigma_rt_t
    d2 = (lnS0k + (mu - v2 / 2.0) * t) / sigma_rt_t
    df = np.exp(-r * t)
    dq = np.exp(-q * t)

    c = s * dq * N(d1) - k * df * N(d2)
    p = k * df * N(-d2) - s * dq * N(-d1)

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and s <= h:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and s >= h:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and s >= h:
        return 0.0
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and s <= h:
        return 0.0
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and s <= h:
        return c
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value and s >= h:
        return c
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value and s >= h:
        return p
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and s <= h:
        return p

    num_observations = 1 + t * nobs

    # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
    # Adjusts the barrier for discrete and not continuous observations
    h_adj = h
    t = t / num_observations

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    else:
        raise FinError("Unknown barrier option type." +
                       str(option_type))

    h = h_adj

    if abs(v) < 1e-5:
        v = 1e-5

    ll = (mu + v*v / 2.0) / v2
    y = np.log(h * h / (s * k)) / sigma_rt_t + ll * sigma_rt_t
    x1 = np.log(s / h) / sigma_rt_t + ll * sigma_rt_t
    y1 = np.log(h / s) / sigma_rt_t + ll * sigma_rt_t
    hOverS = h / s

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        if h >= k:
            c_do = s * dq * N(x1) - k * df * N(x1 - sigma_rt_t) \
                   - s * dq * pow(hOverS, 2.0 * ll) * N(y1) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * N(y1 - sigma_rt_t)
            price = c_do
        else:
            c_di = s * dq * pow(hOverS, 2.0 * ll) * N(y) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * N(y - sigma_rt_t)
            price = c - c_di
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        if h <= k:
            c_di = s * dq * pow(hOverS, 2.0 * ll) * N(y) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * N(y - sigma_rt_t)
            price = c_di
        else:
            c_do = s * dq * N(x1) \
                   - k * df * N(x1 - sigma_rt_t) \
                   - s * dq * pow(hOverS, 2.0 * ll) * N(y1) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * N(y1 - sigma_rt_t)
            price = c - c_do
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        if h >= k:
            c_ui = s * dq * N(x1) - k * df * N(x1 - sigma_rt_t) \
                   - s * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   (N(-y + sigma_rt_t) - N(-y1 + sigma_rt_t))
            price = c_ui
        else:
            price = c
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        if h > k:
            c_ui = s * dq * N(x1) - k * df * N(x1 - sigma_rt_t) \
                   - s * dq * pow(hOverS, 2.0 * ll) * (N(-y) - N(-y1)) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   (N(-y + sigma_rt_t) - N(-y1 + sigma_rt_t))
            price = c - c_ui
        else:
            price = 0.0
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        if h > k:
            p_ui = -s * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigma_rt_t)
            price = p_ui
        else:
            p_uo = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_rt_t) \
                   + s * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   N(-y1 + sigma_rt_t)
            price = p - p_uo
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        if h >= k:
            p_ui = -s * dq * pow(hOverS, 2.0 * ll) * N(-y) \
                   + k * df * pow(hOverS, 2.0 * ll - 2.0) * N(-y + sigma_rt_t)
            price = p - p_ui
        else:
            p_uo = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_rt_t) \
                   + s * dq * pow(hOverS, 2.0 * ll) * N(-y1) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   N(-y1 + sigma_rt_t)
            price = p_uo
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        if h >= k:
            price = 0.0
        else:
            p_di = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_rt_t) \
                   + s * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   (N(y - sigma_rt_t) - N(y1 - sigma_rt_t))
            price = p - p_di
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        if h >= k:
            price = p
        else:
            p_di = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_rt_t) \
                   + s * dq * pow(hOverS, 2.0 * ll) * (N(y) - N(y1)) \
                   - k * df * pow(hOverS, 2.0 * ll - 2.0) * \
                   (N(y - sigma_rt_t) - N(y1 - sigma_rt_t))
            price = p_di
    else:
        raise FinError("Unknown barrier option type." +
                       str(option_type))

    return price


###############################################################################
