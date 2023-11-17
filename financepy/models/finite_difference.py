from copy import deepcopy
from functools import partial

import numpy as np

from ..utils.math import band_matrix_multiplication
from ..utils.math import solve_tridiagonal_matrix
from ..utils.math import transpose_tridiagonal_matrix
from financepy.utils.global_types import OptionTypes


# @njit
def dx(x, wind=0):
    # Intermediate rows
    # Note: As first and last rows are handled separately
    # (at the end of this method)
    # we can use numpy roll without worrying about the end values
    dxl = (x - np.roll(x, 1))
    dxu = (np.roll(x, -1) - x)
    if wind < 0:
        # (-1/dxl, 1/dxl, 0)
        out = np.array((-1 / dxl, 1 / dxl, np.zeros_like(dxl))).T
    elif wind == 0:
        intermediate_rows = np.array(
            (- dxu / dxl,
             dxu / dxl - dxl / dxu,
             dxl / dxu
             )
        ) / (dxl + dxu)
        out = intermediate_rows.T
    else:
        # (0, -1/dxu, 1/dxu)
        out = np.array((np.zeros_like(dxu), -1 / dxu, 1 / dxu)).T

    # First row
    if wind >= 0:
        out[0] = (0, -1, 1)
        out[0] /= (x[1] - x[0])
    else:
        out[0] = (0, 0, 0)

    # Last row
    if wind <= 0:
        out[-1] = (-1, 1, 0)
        out[-1] /= (x[-1] - x[-2])
    else:
        out[-1] = (0, 0, 0)

    return out


def dxx(x):
    # Intermediate rows
    # Note: As first and last rows are handled separately
    # (they overwritten at end of this method),
    # we can use numpy roll without worrying about the end values
    dxl = (x - np.roll(x, 1))
    dxu = (np.roll(x, -1) - x)
    intermediate_rows = np.array(
        [2 / dxl, -(2 / dxl + 2 / dxu), 2 / dxu]) / (dxu + dxl)
    out = intermediate_rows.T

    # First row
    out[0] = (0, 0, 0)
    # Last row
    out[-1] = (0, 0, 0)

    return out


def calculate_fd_matrix(x, r, mu, var, dt, theta, wind=0):
    """
    1d finite difference solution for pdes of the form

    0 = dV/dt + A V

    A = -risk_free_rate + mu d/dx + 1/2 var d2/dx2

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


def fd_roll_backwards(res, theta, Ai=None, Ae=None):
    # TODO Test for more than one vector
    Ai = np.array([]) if Ai is None else Ai
    Ae = np.array([]) if Ae is None else Ae
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


def fd_roll_forwards(res, theta, Ai=None, Ae=None):
    Ai = np.array([]) if Ai is None else Ai
    Ae = np.array([]) if Ae is None else Ae
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

###############################################################################


def option_payoff(s, strike, smooth, dig, option_type):

    if isinstance(option_type, OptionTypes):
        option_type = option_type.value

    # Generate middle values (i.e. not first or last, which are
    # overwritten later)
    if not smooth:
        if dig:
            res = digital(s, strike)
        else:
            res = s - strike
            # Set negative values to zero
            res[res < 0] = 0
    else:
        # Set lower and upper bound for s.
        # Note: As we are only interested in elements 1:-1 here,
        # we can use roll without worrying about end values.
        sl = 0.5 * (np.roll(s, 1) + s)
        su = 0.5 * (s + np.roll(s, -1))

        # Define the curve, fix the strike_price value, and map into res
        func = smooth_digital if dig else smooth_call
        func = partial(func, strike=strike)
        res = list(map(func, sl, su))

    # Handle first and last values separately
    res[0] = digital(s[0], strike) if dig else max(0, s[0] - strike)
    res[-1] = digital(s[-1], strike) if dig else max(0, s[-1] - strike)

    # Invert for put options
    if option_type in {OptionTypes.AMERICAN_PUT.value,
                       OptionTypes.EUROPEAN_PUT.value}:
        res = 1 - res if dig else res - (s - strike)

    return np.atleast_2d(res)


###############################################################################

def black_scholes_fd(spot_price, volatility, time_to_expiry,
                     strike_price, risk_free_rate,
                     dividend_yield, option_type,
                     num_time_steps=None, num_samples=2000,
                     num_std=5, theta=0.5, wind=0,
                     digital=False,
                     smooth=False, update=False):
    if isinstance(option_type, OptionTypes):
        option_type = option_type.value

    # Define grid
    std = volatility * (time_to_expiry ** 0.5)
    xu = num_std * std
    xl = -xu
    d_x = (xu - xl) / max(1, num_samples)
    num_samples = 1 if num_samples <= 0 or xl == xu else num_samples + 1

    # Calculate the drift
    mu = risk_free_rate - dividend_yield

    # Create sample set s
    s = spot_price * np.exp(xl + d_x * np.arange(0, num_samples))

    # Generate the option payoff to be fitted
    payoff = option_payoff(s, strike_price, smooth, digital, option_type)

    # time steps
    num_steps = num_time_steps or num_samples // 2
    dt = time_to_expiry / max(1, num_steps)

    # Make time series for interest rate, drift, and variance
    r_ = np.zeros(num_samples) + risk_free_rate
    mu_ = mu * s
    var_ = (s * volatility) ** 2

    # Initialise implicit and explicit matricies
    Ai = np.array([])
    Ae = np.array([])

    # Store original res as res0
    res = deepcopy(payoff)

    for h in range(num_steps):
        if update or h == 0:
            # Explicit case
            if theta != 1:
                Ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1-theta, wind)
            # Implicit case
            if theta != 0:
                Ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

        res = fd_roll_backwards(res, theta, Ai=Ai, Ae=Ae)

        if option_type in {OptionTypes.AMERICAN_CALL.value,
                           OptionTypes.AMERICAN_PUT.value}:
            idx = res[0] < payoff[0]
            res[0][idx] = payoff[0][idx]

    return res[0][num_samples // 2]
