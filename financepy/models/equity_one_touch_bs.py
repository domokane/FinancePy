##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit
import numba as nb

########################################################################################
# TODO: Implement Sobol random numbers
# TODO: Improve convergence
########################################################################################

DEBUG_MODE = False


@njit(fastmath=True, cache=True)
def barrier_pay_one_at_hit_pv_down(s, hh, r, dt):
    """Pay $1 if the stock crosses the barrier hh from above. PV payment."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] <= hh:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def barrier_pay_one_at_hit_pv_up(s, hh, r, dt):
    """Pay $1 if the stock crosses the barrier hh from below. PV payment."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 0

        for it in range(0, num_time_steps):
            if s[ip][it] >= hh:
                hit_time = dt * it
                v = np.exp(-r * hit_time)
                hit_flag = 1
                break

        pv = pv + v * hit_flag

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def barrier_pay_asset_at_expiry_down_out(s, hh):
    """Pay S(T) at expiry if the down barrier has not been touched."""
    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] <= hh:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


########################################################################################


@njit(fastmath=True, cache=True)
def barrier_pay_asset_at_expiry_up_out(s, hh):
    """Pay S(T) at expiry if the up barrier has not been touched."""

    num_paths, num_time_steps = s.shape
    pv = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 1

        for it in range(0, num_time_steps):
            if s[ip][it] >= hh:
                hit_flag = 0
                break

        pv = pv + hit_flag * s[ip][num_time_steps - 1]

    pv = pv / num_paths
    return pv


########################################################################################

@njit(fastmath=True, cache=True)
def barrier_pay_asset_at_expiry_down(s, hh):
    """Pay S(T) at expiry if a down barrier has been touched."""

    num_paths, num_time_steps = s.shape
    payoff = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 0

        for it in range(num_time_steps):
            if s[ip, it] <= hh:
                hit_flag = 1
                break

        payoff += hit_flag * s[ip, num_time_steps - 1]

    payoff = payoff / num_paths
    return payoff


########################################################################################


@njit(fastmath=True, cache=True)
def barrier_pay_asset_at_expiry_up(s, hh):
    """Pay S(T) at expiry if an up barrier has been touched."""

    num_paths, num_time_steps = s.shape
    payoff = 0.0

    for ip in nb.prange(num_paths):
        hit_flag = 0

        for it in range(num_time_steps):
            if s[ip, it] >= hh:
                hit_flag = 1
                break

        payoff += hit_flag * s[ip, num_time_steps - 1]

    payoff = payoff / num_paths
    return payoff
