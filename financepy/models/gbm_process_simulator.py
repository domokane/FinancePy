##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit  # , float64, int64
from ..utils.math import cholesky
from ..utils.error import FinError

###############################################################################
# WE SIMULATE GEOMETRIC BROWNIAN MOTION WITH DIFFERENT CUTS ACROSS 3D SPACE
# OF ASSETS, PATHS AND TIME STEPS.
###############################################################################

# Seed the random generator globally for reproducibility
np.random.seed(42)


@njit
def get_paths_times(
    num_paths, num_time_steps, t, mu, stock_price, volatility, seed
):
    """Get the simulated GBM process for a single asset with even num paths and
    time steps. Inputs include the number of time steps, paths, the drift mu,
    stock price, volatility and a seed.

    Parameters:
    - num_paths: int, number of paths (must be even)
    - num_time_steps: int, number of time steps
    - t: float, total time in years
    - mu: float, annual drift (mean return)
    - stock_price: float, initial stock price
    - volatility: float, annual volatility
    - seed: int, random seed for reproducibility
    """

    if num_paths % 2 == 0:
        num_paths_even = num_paths // 2
    else:
        raise FinError("Number of paths must be an even number.")

    np.random.seed(seed)
    dt = t / num_time_steps
    vsqrt_dt = volatility * np.sqrt(dt)
    m = np.exp((mu - (volatility**2) / 2) * dt)

    t_all = np.linspace(0.0, t, num_time_steps + 1)
    s_all = np.empty((num_paths, num_time_steps + 1))
    s_all[:, 0] = stock_price

    for it in range(1, num_time_steps + 1):
        g_1d = np.random.standard_normal((num_paths_even))
        for ip in range(0, num_paths_even):
            w = np.exp(g_1d[ip] * vsqrt_dt)
            ip_start = ip * 2
            s_all[ip_start, it] = s_all[ip_start, it - 1] * m * w
            s_all[ip_start + 1, it] = s_all[ip_start + 1, it - 1] * m / w

    return t_all, s_all


###############################################################################


@njit
def get_assets_paths_times(
    num_assets,
    num_paths,
    num_time_steps,
    t,
    mus,
    stock_prices,
    volatilities,
    corr_matrix,
    seed,
):
    """Get the simulated GBM process for a number of assets and paths and num
    time steps. Inputs include the number of assets, paths, the vector of mus,
    stock prices, volatilities, a correlation matrix and a seed.

    Parameters:
    - num_assets: int, number of assets
    - num_paths: int, number of paths (must be even)
    - num_time_steps: int, number of time steps
    - t: float, total time in years
    - mus: ndarray, annual drift rates for each asset
    - stock_prices: ndarray, initial stock prices for each asset
    - volatilities: ndarray, annual volatilities for each asset
    - corr_matrix: ndarray, correlation matrix between assets
    - seed: int, random seed for reproducibility
    """

    if num_paths % 2 == 0:
        num_paths_even = num_paths // 2
    else:
        raise FinError("Number of paths must be an even number.")

    if stock_prices.shape[0] != num_assets:
        raise FinError("Stock price vector incorrect size.")

    if volatilities.shape[0] != num_assets:
        raise FinError("Volatilities vector incorrect size.")

    if mus.shape[0] != num_assets:
        raise FinError("Drift mu vector incorrect size.")

    if (
        corr_matrix.shape[0] != num_assets
        and corr_matrix.shape[1] != num_assets
    ):
        raise FinError("Correlation matrix incorrect size.")

    np.random.seed(seed)
    dt = t / num_time_steps
    vsqrt_dts = volatilities * np.sqrt(dt)
    m = np.exp((mus - volatilities * volatilities / 2.0) * dt)

    s_all = np.empty((num_assets, num_paths, num_time_steps + 1))
    t_all = np.linspace(0, t, num_time_steps + 1)

    g = np.random.standard_normal(
        (num_paths_even, num_time_steps + 1, num_assets)
    )
    c = cholesky(corr_matrix)
    g_corr = np.empty((num_paths_even, num_time_steps + 1, num_assets))

    # or use g_corr = np.einsum('ijk,kl->ijl', g, c)

    # Calculate the Cholesky dot product
    for ip in range(0, num_paths_even):
        for it in range(0, num_time_steps + 1):
            for ia in range(0, num_assets):
                g_corr[ip][it][ia] = 0.0
                for ib in range(0, num_assets):
                    g_corr[ip][it][ia] += g[ip][it][ib] * c[ia][ib]

    for ia in range(0, num_assets):
        s_all[ia, :, 0] = stock_prices[ia]

    for ip in range(0, num_paths_even):
        ip_start = ip * 2
        for it in range(1, num_time_steps + 1):
            for ia in range(0, num_assets):
                z = g_corr[ip, it, ia]
                w = np.exp(z * vsqrt_dts[ia])
                v = m[ia]
                s_all[ia, ip_start, it] = s_all[ia, ip_start, it - 1] * v * w
                s_all[ia, ip_start + 1, it] = (
                    s_all[ia, ip_start + 1, it - 1] * v / w
                )

    return t_all, s_all


###############################################################################


@njit
def get_assets_paths(
    num_assets,
    num_paths,
    t,
    mus,
    stock_prices,
    volatilities,
    corr_matrix,
    seed,
):
    """Get the simulated GBM process for a number of assets and paths for one
    time step. Inputs include the number of assets, paths, the vector of mus,
    stock prices, volatilities, a correlation matrix and a seed.

    Parameters:
    - num_assets: int, number of assets
    - num_paths: int, number of paths (must be even)
    - t: float, total time in years
    - mus: ndarray, annual drift rates for each asset
    - stock_prices: ndarray, initial stock prices for each asset
    - volatilities: ndarray, annual volatilities for each asset
    - corr_matrix: ndarray, correlation matrix between assets
    - seed: int, random seed for reproducibility
    """

    if num_paths % 2 == 0:
        num_paths_even = num_paths // 2
    else:
        raise FinError("Number of paths must be an even number.")

    if stock_prices.shape[0] != num_assets:
        raise FinError("Stock price vector incorrect size.")

    if volatilities.shape[0] != num_assets:
        raise FinError("Volatilities vector incorrect size.")

    if mus.shape[0] != num_assets:
        raise FinError("Drift mu vector incorrect size.")

    if (
        corr_matrix.shape[0] != num_assets
        or corr_matrix.shape[1] != num_assets
    ):
        raise FinError("Correlation matrix incorrect size.")

    np.random.seed(seed)
    vsqrt_dts = volatilities * np.sqrt(t)
    m = np.exp((mus - volatilities * volatilities / 2.0) * t)
    s_all = np.empty((num_assets, 2 * num_paths_even))

    num_time_steps = 1
    t_all = np.linspace(0, t, num_time_steps + 1)

    g = np.random.standard_normal((num_paths_even, num_assets))
    c = cholesky(corr_matrix)
    g_corr = np.empty((num_paths_even, num_assets))

    # Calculate the dot product
    for ia in range(0, num_assets):
        for ip in range(0, num_paths_even):
            g_corr[ip][ia] = 0.0
            for ib in range(0, num_assets):
                g_corr[ip][ia] += g[ip][ib] * c[ia][ib]

    for ia in range(0, num_assets):
        for ip in range(0, num_paths_even):
            z = g_corr[ip, ia]
            w = np.exp(z * vsqrt_dts[ia])
            ip_start = ip * 2
            s_all[ia, ip_start] = stock_prices[ia] * m[ia] * w
            s_all[ia, ip_start + 1] = stock_prices[ia] * m[ia] / w

    return t_all, s_all


###############################################################################
