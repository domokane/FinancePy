from enum import Enum
from copy import copy
from functools import partial

import numpy as np

from ..utils.math import band_matrix_multiplication, solve_tridiagonal_matrix, transpose_tridiagonal_matrix
from ..utils.global_vars import gDaysInYear


class PUT_CALL(Enum):
    PUT = -1
    CALL = 1


class AMER_EURO(Enum):
    EURO = 0
    AMER = 1


def dx(x, wind=0):
    n = len(x) - 1
    out = np.zeros((n + 1, 3))

    dxu = x[1] - x[0]

    # First row
    if wind >= 0:
        out[0] = (-1, 0, 1)
        out[0] /= dxu

    # Intermediate rows
    # Note: As first and last rows are handled separately,
    # we can use numpy roll without worrying about the end values
    dxl = (x - np.roll(x, 1))[1:-1]
    dxu = (np.roll(x, -1) - x)[1:-1]
    if wind < 0:
        # (-1/dxl, 1/dxl, 0)
        out[1:-1] = np.array((-1 / dxl, 1 / dxl, np.zeros_like(dxl))).T
    elif wind == 0:
        intermediate_rows = np.array(
            (- dxu / dxl,
             dxu / dxl - dxl / dxu,
             dxl / dxu
             )
        ) / (dxl + dxu)
        out[1:-1] = intermediate_rows.T
    else:
        # (0, -1/dxu, 1/dxu)
        out[1:-1] = np.array((np.zeros_like(dxu), -1 / dxu, 1 / dxu)).T

    # Last row
    if wind <= 0:
        dxl = x[n] - x[n - 1]
        out[n] = (-1, 1, 0)
        out[n] /= dxl
    else:
        out[n] = (0, 0, 0)

    return out


def dxx(x):
    n = len(x) - 1
    out = np.zeros((n + 1, 3))

    # First row
    out[0] = (0, 0, 0)

    # Intermediate rows
    # Note: As first and last rows are handled separately,
    # we can use numpy roll without worrying about the end values
    dxl = (x - np.roll(x, 1))[1:-1]
    dxu = (np.roll(x, -1) - x)[1:-1]
    intermediate_rows = np.array(
        [2 / dxl,
        -(2 / dxl + 2 / dxu),
        2 / dxu]
    ) / (dxu + dxl)
    out[1:-1] = intermediate_rows.T

    # Last row
    out[n] = (0, 0, 0)

    return out


def calculate_fd_matrix(x, r, mu, var, dt, theta, wind=0):
    """
    1d finite difference solution for pdes of the form

    0 = dV/dt + A V

    A = -risk_free_rate + mu d/dx + 1/2 var d^2/dx^2

    using the theta scheme

    [1-theta dt A] V(t) = [1 + (1-theta) dt A] V(t+dt)
    """
    if dt == 0:
        raise ValueError("Timestep length dt must be non-zero")
    if theta == 0:
        raise ValueError("Theta must be non-zero")

    # Calculate the finite differences for the first and second derivatives
    Dxx = dxx(x)

    if wind == 0:
        Dx = dx(x, 0)
    elif wind < 0:
        Dxd = dx(x, -1)
        Dx = Dxd
    elif wind == 1:
        Dxu = dx(x, 1)
        Dx = Dxu
    elif wind > 1:
        # use Dxd when mu < 0 and Dxu otherwise
        Dxd = dx(x, -1)
        Dxu = dx(x, 1)
        Dx = np.zeros((len(x), 1)) + Dxu
        Dx[mu[0] < 0] = Dxd

    # Ensure mu and var have correct dimensions
    mu = np.atleast_2d(mu)
    var = np.atleast_2d(var)

    # Calculate matrix
    mm = Dx.shape[1] // 2  # integer division
    A = dt * theta * (mu.T * Dx + 0.5 * var.T * Dxx)
    A[:, mm] += 1 - dt * theta * r

    return A


def fd_roll_backwards(res, theta, Ai=np.array([]), Ae=np.array([])):
    num_vectors = len(res)
    mm = 1

    # Explicit case
    if theta != 1:
        for k in range(num_vectors):
            res[k] = band_matrix_multiplication(Ae, mm, mm, res[k])

    # Implicit case
    if theta != 0:
        for k in range(num_vectors):
            res[k] = solve_tridiagonal_matrix(Ai, res[k])

    return res


def fd_roll_forwards(res, theta, Ai=np.array([]), Ae=np.array([])):
    num_vectors = len(res)
    mm = num_vectors // 2  # integer division

    # Implicit case
    if theta != 0:
        Ai = transpose_tridiagonal_matrix(Ai)
        for k in range(num_vectors):
            res[k] = solve_tridiagonal_matrix(Ai, res[k])

    # Explicit case
    if theta != 1:
        Ae = transpose_tridiagonal_matrix(Ae)
        for k in range(num_vectors):
            res[k] = band_matrix_multiplication(Ae, mm, mm, res[k])

    return res


def smooth_digital(xl, xu, strike):
    if xu <= strike:
        return 0
    elif strike <= xl:
        return 1
    else:
        return (xu - strike) / (xu - xl)


def digital(x, strike):
    return 0.5 * (np.sign(x - strike) + 1)


def smooth_call(xl, xu, strike):
    if xu <= strike:
        return 0
    elif strike <= xl:
        return 0.5 * (xl + xu) - strike
    else:
        return 0.5 * (xu - strike) ** 2 / (xu - xl)


def initial_curve(s, strike, smooth, dig, pc):
    num_samples = len(s)
    res = np.zeros(num_samples)

    # Handle first and last values separately
    res[0] = digital(s[0], strike) if dig else max(0, s[0] - strike)
    res[-1] = digital(s[-1], strike) if dig else max(0, s[-1] - strike)

    # Generate middle values (i.e. not first or last)
    if not smooth:
        if dig:
            res[1:-1] = digital(s[1:-1], strike)
        else:
            res[1:-1] = s[1:-1] - strike
            # Set negative values to zero
            res[1:-1][res[1:-1] < 0] = 0
    else:
        # Set lower and upper bound for s.
        # Note: As we are only interested in elements 1:-1 here,
        # we can use roll without worrying about end values.
        sl = 0.5 * (np.roll(s, 1) + s)
        su = 0.5 * (s + np.roll(s, -1))

        # Define the curve, fix the strike_price value, and map into res
        func = smooth_digital if dig else smooth_call
        func = partial(func, strike=strike)
        res[1:-1] = list(map(func, sl[1:-1], su[1:-1]))

    # Invert for put options
    if pc == PUT_CALL.PUT.value:
        res = 1 - res if dig else res - (s - strike)

    return np.atleast_2d(res)


def black_scholes_finite_difference(stock_price, risk_free_rate, mu, sigma, expiry_date, valuation_date, strike_price,
                                    dig, pc, exercise, smooth, theta, wind, num_std, num_steps, num_samples, update, num_pr):
    time_to_expiry = (expiry_date - valuation_date) / gDaysInYear
    std = sigma * (time_to_expiry ** 0.5)
    xl = -num_std * std
    xu = num_std * std
    d_x = (xu - xl) / max(1, num_samples)
    num_samples = 1 if num_samples <= 0 or xl == xu else num_samples + 1

    # Create sample set s
    s = np.zeros(num_samples)
    s[0] = stock_price * np.exp(xl)
    ds = np.exp(d_x)
    for i in range(1, num_samples):
        s[i] = s[i - 1] * ds

    # Define the initial curve which will be fitted with each iteration
    res = initial_curve(s, strike_price, smooth, dig, pc)

    # time steps
    dt = time_to_expiry / max(1, num_steps)

    # Make time series
    r_ = np.zeros(num_samples) + risk_free_rate
    mu_ = mu * s
    var_ = (s * sigma) ** 2

    # Store original res if option is American
    if exercise == AMER_EURO.AMER.value:
        res0 = copy(res)

    # repeat
    nump = max(1, num_pr)
    for p in range(nump):
        Ai = np.array([])
        Ae = np.array([])
        for h in range(num_steps):
            if update or h == 0:
                # Explicit case
                if theta != 1:
                    Ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1-theta, wind)
                # Implicit case
                if theta != 0:
                    Ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

            res = fd_roll_backwards(res, theta, Ai=Ai, Ae=Ae)
            if exercise == AMER_EURO.AMER.value:
                idx = res[0] < res0[0]
                res[0][idx] = res0[0][idx]

    return res[0], res[0][num_samples // 2]
