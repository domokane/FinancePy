from enum import Enum
from copy import copy

import numpy as np

from ..utils.math import band_matrix_multiplication, solve_tridiagonal_matrix, transpose_tridiagonal_matrix, sign

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
    for i in range(1, n):
        dxl = x[i] - x[i - 1]
        dxu = x[i + 1] - x[i]
        if wind < 0:
            out[i] = (-1, 1, 0)
            out[i] /= dxl
        elif wind == 0:
            out[i] = (- dxu / dxl,
                      dxu / dxl - dxl / dxu,
                      dxl / dxu)
            out[i] /= (dxl + dxu)
        else:
            out[i] = (0, -1, 1)
            out[i] /= dxu

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
    for i in range(1, n):
        dxl = x[i] - x[i - 1]
        dxu = x[i + 1] - x[i]

        out[i] = (1 / dxl,
                  -(1 / dxl + 1 / dxu),
                  1 / dxu
                  )
        out[i] *= 2.0 / (dxl + dxu)

    # Last row
    out[n] = (0, 0, 0)

    return out


def calculate_fd_matrix(x, r, mu, var, dt, theta, wind=0):
    """
    1d finite difference solution for pdes of the form

    0 = dV/dt + A V

    A = -r + mu d/dx + 1/2 var d^2/dx^2

    using the theta scheme

    [1-theta dt A] V(t) = [1 + (1-theta) dt A] V(t+dt)
    """
    if dt == 0:
        raise ValueError("Timestep length dt must be non-zero")
    if theta == 0:
        raise ValueError("Theta must be non-zero")

    Dxd = dx(x, -1)
    Dxu = dx(x, 1)
    Dxx = dxx(x)
    Dx = dx(x, 0)

    if wind < 0:
        Dx = Dxd
    elif wind == 1:
        Dx = Dxu

    n = len(x)
    m = Dx.shape[1]
    A = np.zeros((n, m))

    mm = m // 2  # integer division
    dtTheta = dt * theta
    for i in range(n):
        if wind > 1:
            Dx = Dxd if mu[i] < 0 else Dxu

        for j in range(m):
            A[i, j] = dtTheta * (mu[i] * Dx[i, j] + 0.5 * var[i] * Dxx[i, j])

        A[i, mm] += 1 - dtTheta * r[i]

    return A


def fd_roll_backwards(x, r, mu, var, dt, res, theta, update, Ai=np.array([]), Ae=np.array([]), wind=0):
    num_vectors = len(res)
    mm = 1

    # Explicit case
    if theta != 1:
        if update:
            Ai = calculate_fd_matrix(x, r, mu, var, dt, 1 - theta, wind)
        for k in range(num_vectors):
            res[k] = band_matrix_multiplication(Ai, mm, mm, res[k])

    # Implicit case
    if theta != 0:
        if update:
            Ae = calculate_fd_matrix(x, r, mu, var, -dt, theta, wind)
        for k in range(num_vectors):
            res[k] = solve_tridiagonal_matrix(Ae, res[k])

    return res, Ai, Ae


def fd_roll_forwards(x, r, mu, var, dt, res, theta, update, Ai=np.array([]), Ae=np.array([]), wind=0):
    num_vectors = len(res)
    mm = num_vectors // 2  # integer division

    # Implicit case
    if theta != 0:
        if update:
            Ai = calculate_fd_matrix(x, r, mu, var, -dt, theta, wind)
            Ai = transpose_tridiagonal_matrix(Ai)
        for k in range(num_vectors):
            res[k] = solve_tridiagonal_matrix(Ai, res[k])

    # Explicit case
    if theta != 1:
        if update:
            Ae = calculate_fd_matrix(x, r, mu, var, dt, 1 - theta, wind)
            Ae = transpose_tridiagonal_matrix(Ae)
        for k in range(num_vectors):
            res[k] = band_matrix_multiplication(Ae, mm, mm, res[k])

    return res, Ai, Ae


def smooth_digital(xl, xu, strike):
    if xu <= strike:
        res = 0
    elif strike <= xl:
        res = 1
    else:
        res = (xu - strike) / (xu - xl)

    return res


def smooth_call(xl, xu, strike):
    if xu <= strike:
        res = 0
    elif strike <= xl:
        res = 0.5 * (xl + xu) - strike
    else:
        res = 0.5 * (xu - strike) ** 2 / (xu - xl)

    return res


def black_scholes_finite_difference(s0, r, mu, sigma, expiry, strike, dig, pc, ea, smooth, theta, wind,
                                    num_std, num_t, num_s, update, num_pr):
    t = max(0, expiry)
    std = sigma * (t ** 0.5)
    xl = -num_std * std
    xu = num_std * std
    d_x = (xu - xl) / max(1, num_s)
    nums = num_s
    if nums <= 0 or xl == xu:
        nums = 1
    else:
        nums += 1

    # Create sample set s
    s = np.zeros(nums)
    s[0] = s0 * np.exp(xl)
    ds = np.exp(d_x)

    for i in range(1, nums):
        s[i] = s[i - 1] * ds
    res = np.zeros(nums)
    for i in range(nums):
        if smooth == 0 or i == 0 or i == nums - 1:
            if dig:
                res[i] = 0.5 * (sign(s[i] - strike) + 1)
            else:
                res[i] = max(0, s[i] - strike)
        else:
            sl = 0.5 * (s[i - 1] + s[i]);
            su = 0.5 * (s[i] + s[i + 1]);
            if dig:
                res[i] = smooth_digital(sl, su, strike)
            else:
                res[i] = smooth_call(sl, su, strike)

        if pc == PUT_CALL.PUT.value:
            if dig:
                res[i] = 1 - res[i]
            else:
                res[i] -= (s[i] - strike)

    # time steps
    numt = max(0, num_t)
    dt = t / max(1, numt)

    # repeat
    res = np.array([res])
    res0 = copy(res)
    nump = max(1, num_pr)
    for p in range(nump):
        r_ = np.zeros(nums) + r
        mu_ = mu * s
        var_ = (s * sigma) ** 2
        Ai = np.array([])
        Ae = np.array([])
        for h in range(numt):
            res, Ai, Ae = fd_roll_backwards(
                s, r_, mu_, var_, dt, res, theta, update=update or h == 0, Ai=Ai, Ae=Ae, wind=wind)
            if ea == AMER_EURO.AMER.value:
                for i in range(nums):
                    res[0][i] = max(res[0][i], res0[0][i])

    return res[0], res[0][nums // 2]
