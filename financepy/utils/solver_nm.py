########################################################################################
# Copyright (C) 2020 Saeed Amen, Dominic O'Kane
########################################################################################

from collections import namedtuple
import numpy as np
from numba import njit

########################################################################################
# Result structure
results = namedtuple("results", "x fun success nit final_simplex")

########################################################################################
# Default tolerances
_XTOL = 2e-12
_RTOL = 4 * np.finfo(float).eps

########################################################################################


@njit(fastmath=True)
def nelder_mead(
    fun,
    x0,
    bounds=np.array([[], []]).T,
    args=(),
    tol_f=1e-10,
    tol_x=1e-10,
    max_iter=1000,
    roh=1.0,
    chi=2.0,
    v=0.5,
    sigma=0.5,
):
    """
    Minimize a scalar-valued function using Nelder-Mead method.

    Parameters
    ----------
    fun : callable
        Function to minimize, must accept (x, *args)
    x0 : ndarray
        Initial guess, shape (n,)
    bounds : ndarray, optional
        Shape (n,2) of (min,max) for each variable
    args : tuple
        Extra arguments to fun
    tol_f, tol_x : float
        Tolerances for convergence
    max_iter : int
        Maximum iterations
    roh, chi, v, sigma : float
        Nelder-Mead parameters

    Returns
    -------
    results namedtuple
    """
    n = x0.size
    vertices = _initialize_simplex(x0, bounds)
    _check_params(roh, chi, v, sigma, bounds, n)

    f_val = np.empty(n + 1, dtype=np.float64)
    for i in range(n + 1):
        f_val[i] = fun(vertices[i], *args)

    sort_ind = f_val.argsort()
    lv_ratio = 1.0
    nit = 0
    x_bar = vertices[sort_ind[:n]].sum(axis=0) / n
    sigma_n = sigma**n
    rohv = roh * v
    rohchi = roh * chi

    while True:
        shrink = False
        best_val_idx = sort_ind[0]
        worst_val_idx = sort_ind[n]

        term_f = f_val[worst_val_idx] - f_val[best_val_idx] < tol_f
        term_x = lv_ratio < tol_x
        fail = nit >= max_iter

        if term_f or term_x or fail:
            break

        # Step 2: Reflection
        x_r = x_bar + roh * (x_bar - vertices[worst_val_idx])
        f_r = fun(x_r, *args)

        if f_val[best_val_idx] <= f_r < f_val[sort_ind[n - 1]]:
            vertices[worst_val_idx] = x_r
            lv_ratio *= roh

        elif f_r < f_val[best_val_idx]:
            x_e = x_bar + chi * (x_r - x_bar)
            f_e = fun(x_e, *args)
            if f_e < f_r:
                vertices[worst_val_idx] = x_e
                lv_ratio *= rohchi
            else:
                vertices[worst_val_idx] = x_r
                lv_ratio *= roh
        else:
            # Contraction
            if f_r < f_val[worst_val_idx]:
                x_c = x_bar + v * (x_r - x_bar)
                lv_ratio_update = rohv
            else:
                x_c = x_bar - v * (x_r - x_bar)
                lv_ratio_update = v

            f_c = fun(x_c, *args)
            if f_c < min(f_r, f_val[worst_val_idx]):
                vertices[worst_val_idx] = x_c
                lv_ratio *= lv_ratio_update
            else:
                # Shrink
                shrink = True
                best_vertex = vertices[best_val_idx].copy()
                for i in sort_ind[1:]:
                    vertices[i] = best_vertex + sigma * (vertices[i] - best_vertex)
                    f_val[i] = fun(vertices[i], *args)
                sort_ind[1:] = f_val[sort_ind[1:]].argsort() + 1
                x_bar = (
                    best_vertex
                    + sigma * (x_bar - best_vertex)
                    + (vertices[worst_val_idx] - vertices[sort_ind[n]]) / n
                )
                lv_ratio *= sigma_n

        if not shrink:
            f_val[worst_val_idx] = fun(vertices[worst_val_idx], *args)
            for i, j in enumerate(sort_ind):
                if f_val[worst_val_idx] < f_val[j]:
                    sort_ind[i + 1 :] = sort_ind[i:-1]
                    sort_ind[i] = worst_val_idx
                    break
            x_bar += (vertices[worst_val_idx] - vertices[sort_ind[n]]) / n

        nit += 1

    best_x = vertices[sort_ind[0]]
    best_f = f_val[sort_ind[0]]
    success = not fail

    return results(
        x=best_x, fun=best_f, success=success, nit=nit, final_simplex=vertices
    )


########################################################################################
# Helper functions


@njit(cache=True, fastmath=True)
def _initialize_simplex(x0, bounds):
    n = x0.size
    vertices = np.empty((n + 1, n), dtype=np.float64)
    vertices[:] = x0
    nonzdelt = 0.05
    zdelt = 0.00025
    for i in range(n):
        if vertices[i + 1, i] != 0.0:
            vertices[i + 1, i] *= 1 + nonzdelt
        else:
            vertices[i + 1, i] = zdelt
    return vertices


@njit(cache=True, fastmath=True)
def _check_params(rho, chi, v, sigma, bounds, n):
    if rho <= 0:
        raise ValueError("rho must be > 0")
    if chi <= max(1, rho):
        raise ValueError("chi must be > max(1, rho)")
    if not 0 < v < 1:
        raise ValueError("v must be between 0 and 1")
    if not 0 < sigma < 1:
        raise ValueError("sigma must be between 0 and 1")
    if not (bounds.shape == (0, 2) or bounds.shape == (n, 2)):
        raise ValueError("bounds shape invalid")
    if (np.atleast_2d(bounds)[:, 0] > np.atleast_2d(bounds)[:, 1]).any():
        raise ValueError("lower bound > upper bound")
