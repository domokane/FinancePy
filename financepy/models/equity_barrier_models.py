# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import float64, int64, vectorize

from ..utils.error import FinError
from ..utils.global_types import EquityBarrierTypes
from ..utils.math import normcdf


# Calculates the Barrier option price using an Analytical Approach
# and the Black Scholes Model

########################################################################################


@vectorize(
    [
        float64(
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            int64,
            int64,
        )
    ],
    fastmath=True,
    cache=True,
)
def value_barrier(t, k, h, s, r, q, v, opt_type, nobs):
    """This values a single barrier option. Because of its structure it cannot
    easily be vectorised which is why it has been wrapped.
    # number of observations per year
    """
    ln_s0_k = np.log(s / k)
    sqrt_t = np.sqrt(t)

    sigma_rt_t = v * sqrt_t
    v2 = v * v
    mu = r - q
    d1 = (ln_s0_k + (mu + v2 / 2.0) * t) / sigma_rt_t
    d2 = (ln_s0_k + (mu - v2 / 2.0) * t) / sigma_rt_t
    df = np.exp(-r * t)
    dq = np.exp(-q * t)

    c = s * dq * normcdf(+d1) - k * df * normcdf(+d2)
    p = k * df * normcdf(-d2) - s * dq * normcdf(-d1)

    if s >= h:

        if opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
            return 0.0
        if opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
            return 0.0
        if opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
            return c
        if opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
            return p

    else:

        if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
            return 0.0
        elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
            return 0.0
        elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
            return c
        elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
            return p

    # if 1 == 0:
    #     if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and s <= h:
    #         return 0.0
    #     elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and s >= h:
    #         return 0.0
    #     elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and s >= h:
    #         return 0.0
    #     elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and s <= h:
    #         return 0.0
    #     elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and s <= h:
    #         return c
    #     elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value and s >= h:
    #         return c
    #     elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value and s >= h:
    #         return p
    #     elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and s <= h:
    #         return p

    num_observations = 1 + t * nobs

    # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
    # Adjusts the barrier for discrete and not continuous observations
    h_adj = h
    t = t / num_observations

    if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        h_adj = h * np.exp(0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        h_adj = h * np.exp(-0.5826 * v * np.sqrt(t))
    else:
        raise FinError("Unknown barrier option type." + str(opt_type))

    h = h_adj

    if abs(v) < 1e-5:
        v = 1e-5

    ll = (mu + v * v / 2.0) / v2
    y = np.log(h * h / (s * k)) / sigma_rt_t + ll * sigma_rt_t
    x1 = np.log(s / h) / sigma_rt_t + ll * sigma_rt_t
    y1 = np.log(h / s) / sigma_rt_t + ll * sigma_rt_t
    h_over_s = h / s

    if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        if h >= k:
            c_do = (
                s * dq * normcdf(x1)
                - k * df * normcdf(x1 - sigma_rt_t)
                - s * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
                + k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(y1 - sigma_rt_t)
            )
            price = c_do
        else:
            c_di = s * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
                h_over_s, 2.0 * ll - 2.0
            ) * normcdf(y - sigma_rt_t)
            price = c - c_di
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        if h <= k:
            c_di = s * dq * pow(h_over_s, 2.0 * ll) * normcdf(y) - k * df * pow(
                h_over_s, 2.0 * ll - 2.0
            ) * normcdf(y - sigma_rt_t)
            price = c_di
        else:
            c_do = (
                s * dq * normcdf(x1)
                - k * df * normcdf(x1 - sigma_rt_t)
                - s * dq * pow(h_over_s, 2.0 * ll) * normcdf(y1)
                + k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(y1 - sigma_rt_t)
            )
            price = c - c_do
    elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        if h >= k:
            c_ui = (
                s * dq * normcdf(x1)
                - k * df * normcdf(x1 - sigma_rt_t)
                - s * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
                + k
                * df
                * pow(h_over_s, 2.0 * ll - 2.0)
                * (normcdf(-y + sigma_rt_t) - normcdf(-y1 + sigma_rt_t))
            )
            price = c_ui
        else:
            price = c
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        if h > k:
            c_ui = (
                s * dq * normcdf(x1)
                - k * df * normcdf(x1 - sigma_rt_t)
                - s * dq * pow(h_over_s, 2.0 * ll) * (normcdf(-y) - normcdf(-y1))
                + k
                * df
                * pow(h_over_s, 2.0 * ll - 2.0)
                * (normcdf(-y + sigma_rt_t) - normcdf(-y1 + sigma_rt_t))
            )
            price = c - c_ui
        else:
            price = 0.0
    elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        if h > k:
            p_ui = -s * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
                h_over_s, 2.0 * ll - 2.0
            ) * normcdf(-y + sigma_rt_t)
            price = p_ui
        else:
            p_uo = (
                -s * dq * normcdf(-x1)
                + k * df * normcdf(-x1 + sigma_rt_t)
                + s * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
                - k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(-y1 + sigma_rt_t)
            )
            price = p - p_uo
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        if h >= k:
            p_ui = -s * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y) + k * df * pow(
                h_over_s, 2.0 * ll - 2.0
            ) * normcdf(-y + sigma_rt_t)
            price = p - p_ui
        else:
            p_uo = (
                -s * dq * normcdf(-x1)
                + k * df * normcdf(-x1 + sigma_rt_t)
                + s * dq * pow(h_over_s, 2.0 * ll) * normcdf(-y1)
                - k * df * pow(h_over_s, 2.0 * ll - 2.0) * normcdf(-y1 + sigma_rt_t)
            )
            price = p_uo
    elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        if h >= k:
            price = 0.0
        else:
            p_di = (
                -s * dq * normcdf(-x1)
                + k * df * normcdf(-x1 + sigma_rt_t)
                + s * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
                - k
                * df
                * pow(h_over_s, 2.0 * ll - 2.0)
                * (normcdf(y - sigma_rt_t) - normcdf(y1 - sigma_rt_t))
            )
            price = p - p_di
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        if h >= k:
            price = p
        else:
            p_di = (
                -s * dq * normcdf(-x1)
                + k * df * normcdf(-x1 + sigma_rt_t)
                + s * dq * pow(h_over_s, 2.0 * ll) * (normcdf(y) - normcdf(y1))
                - k
                * df
                * pow(h_over_s, 2.0 * ll - 2.0)
                * (normcdf(y - sigma_rt_t) - normcdf(y1 - sigma_rt_t))
            )
            price = p_di
    else:
        raise FinError("Unknown barrier option type." + str(opt_type))

    return price
