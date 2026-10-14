import numpy as np

from numba import njit

from ..utils.math import band_matrix_multiplication
from ..utils.math import solve_tridiagonal_matrix
from ..utils.math import transpose_tridiagonal_matrix
from ..utils.global_types import OptionTypes

########################################################################################


@njit(cache=True, fastmath=True)
def validate_black_scholes_fd_inputs(
    spot_price,
    volatility,
    time_to_expiry,
    strike_price,
    risk_free_rate,
    dividend_yield,
    opt_type,
    num_steps_per_year,
    num_samples,
    num_std,
    theta,
    wind,
):
    if spot_price <= 0.0:
        raise ValueError("spot_price must be positive")

    if strike_price <= 0.0:
        raise ValueError("strike_price must be positive")

    if risk_free_rate < 0.0:
        raise ValueError("risk_free_rate must be positive")

    if dividend_yield < 0.0:
        raise ValueError("dividend_yield must be positive")

    if volatility < 0.0:
        raise ValueError("volatility must be non-negative")

    if time_to_expiry < 0.0:
        raise ValueError("time_to_expiry must be non-negative")

    if num_samples < 3:
        raise ValueError("num_samples must be at least 3")

    if num_std <= 0.0:
        raise ValueError("num_std must be positive")

    if not (0.0 <= theta <= 1.0):
        raise ValueError("theta must be between 0 and 1")

    if wind not in {-1, 0, 1, 2}:
        raise ValueError("wind must be one of {-1, 0, 1, 2}")

    if num_steps_per_year is not None and num_steps_per_year <= 0:
        raise ValueError("num_steps_per_year must be positive or None")

    valid_types = {
        OptionTypes.EUROPEAN_CALL.value,
        OptionTypes.EUROPEAN_PUT.value,
        OptionTypes.AMERICAN_CALL.value,
        OptionTypes.AMERICAN_PUT.value,
    }

    if opt_type not in valid_types:
        raise ValueError(f"Unsupported option type: {opt_type}")


@njit(cache=True, fastmath=True)
def fn_dx(x, wind=0):

    # Intermediate rows
    # Note: first and last rows are handled separately
    # at the end of this method.
    dxl = x - np.roll(x, 1)
    dxu = np.roll(x, -1) - x

    n = len(x)
    out = np.empty((n, 3), dtype=x.dtype)

    if wind < 0:
        # (-1/dxl, 1/dxl, 0)
        out[:, 0] = -1.0 / dxl
        out[:, 1] = 1.0 / dxl
        out[:, 2] = 0.0

    elif wind == 0:
        denom = dxl + dxu

        out[:, 0] = (-dxu / dxl) / denom
        out[:, 1] = (dxu / dxl - dxl / dxu) / denom
        out[:, 2] = (dxl / dxu) / denom

    else:
        # (0, -1/dxu, 1/dxu)
        out[:, 0] = 0.0
        out[:, 1] = -1.0 / dxu
        out[:, 2] = 1.0 / dxu

    # First row
    if wind >= 0:
        h = x[1] - x[0]
        out[0, 0] = 0.0
        out[0, 1] = -1.0 / h
        out[0, 2] = 1.0 / h
    else:
        out[0, :] = 0.0

    # Last row
    if wind <= 0:
        h = x[-1] - x[-2]
        out[-1, 0] = -1.0 / h
        out[-1, 1] = 1.0 / h
        out[-1, 2] = 0.0
    else:
        out[-1, :] = 0.0

    return out

########################################################################################


@njit(cache=True, fastmath=True)
def fn_dxx(x):

    # Intermediate rows
    # Note: first and last rows are handled separately.
    dxl = x - np.roll(x, 1)
    dxu = np.roll(x, -1) - x

    n = len(x)
    out = np.empty((n, 3), dtype=x.dtype)

    denom = dxu + dxl

    out[:, 0] = (2.0 / dxl) / denom
    out[:, 1] = -(2.0 / dxl + 2.0 / dxu) / denom
    out[:, 2] = (2.0 / dxu) / denom

    # Boundary rows
    out[0, :] = 0.0
    out[-1, :] = 0.0

    return out


########################################################################################


@njit(cache=True, fastmath=True)
def calculate_fd_matrix(x, r, mu, var, dt, theta, wind=0):
    """
    1d finite difference solution for PDEs of the form

    0 = dV/dt + a V

    a = -risk_free_rate + mu d/dx + 1/2 var d2/dx2

    using the theta scheme

    [1-theta dt a] V(t) = [1 + (1-theta) dt a] V(t+dt)
    """

    if dt == 0.0:
        raise ValueError("Timestep length dt must be non-zero")

    if theta == 0.0:
        raise ValueError("Theta must be non-zero")

    dxx = fn_dxx(x)

    if wind == 0:
        dx = fn_dx(x, 0)

    elif wind < 0:
        dx = fn_dx(x, -1)

    elif wind == 1:
        dx = fn_dx(x, 1)

    else:
        # Adaptive upwinding:
        # use downward difference when local drift is negative,
        # upward difference otherwise.

        dxd = fn_dx(x, -1)
        dxu = fn_dx(x, 1)

        dx = dxu.copy()

        for i in range(len(x)):
            if mu[i] < 0.0:
                dx[i, 0] = dxd[i, 0]
                dx[i, 1] = dxd[i, 1]
                dx[i, 2] = dxd[i, 2]

    mm = dx.shape[1] // 2

    # Avoid np.atleast_2d here; explicitly construct the matrix.
    a = np.empty_like(dx)

    for i in range(len(x)):
        a[i, 0] = dt * theta * (
            mu[i] * dx[i, 0]
            + 0.5 * var[i] * dxx[i, 0]
        )

        a[i, 1] = dt * theta * (
            mu[i] * dx[i, 1]
            + 0.5 * var[i] * dxx[i, 1]
        )

        a[i, 2] = dt * theta * (
            mu[i] * dx[i, 2]
            + 0.5 * var[i] * dxx[i, 2]
        )

        a[i, mm] += 1.0 - dt * theta * r[i]

    return a


########################################################################################


@njit(cache=True, fastmath=True)
def fd_roll_backwards(res, theta, ai, ae):

    num_vectors = res.shape[0]
    num_points = res.shape[1]
    mm = 1

    if theta != 1.0:
        for k in range(num_vectors):
            tmp = band_matrix_multiplication(
                ae, mm, mm, res[k]
            )

            for j in range(num_points):
                res[k, j] = tmp[j]

    if theta != 0.0:
        for k in range(num_vectors):
            tmp = solve_tridiagonal_matrix(
                ai, res[k]
            )

            for j in range(num_points):
                res[k, j] = tmp[j]

    return res

########################################################################################


@njit(cache=True, fastmath=True)
def fd_roll_forwards(res, theta, ai=None, ae=None):

    ai = np.array([]) if ai is None else ai
    ae = np.array([]) if ae is None else ae
    num_vectors = len(res)
    mm = num_vectors // 2  # integer division

    # Implicit case
    if theta != 0:
        ai = transpose_tridiagonal_matrix(ai)
        for k in range(num_vectors):
            res[k] = solve_tridiagonal_matrix(ai, res[k])

    # Explicit case
    if theta != 1:
        ae = transpose_tridiagonal_matrix(ae)
        for k in range(num_vectors):
            res[k] = band_matrix_multiplication(ae, mm, mm, res[k])

    return res


########################################################################################


@njit(cache=True, fastmath=True)
def smooth_digital(xl, xu, strike):

    if xu <= strike:
        return 0
    elif strike <= xl:
        return 1

    return (xu - strike) / (xu - xl)


########################################################################################


@njit(cache=True, fastmath=True)
def digital(x, strike):

    return 0.5 * (np.sign(x - strike) + 1)


########################################################################################


@njit(cache=True, fastmath=True)
def smooth_call(xl, xu, strike):

    if xu <= strike:
        return 0
    elif strike <= xl:
        return 0.5 * (xl + xu) - strike

    return 0.5 * (xu - strike) ** 2 / (xu - xl)


########################################################################################

@njit(cache=True, fastmath=True)
def option_payoff(s, strike, smooth, dig: bool, opt_type_int):

    is_put = opt_type_int in {
        OptionTypes.AMERICAN_PUT.value,
        OptionTypes.EUROPEAN_PUT.value,
    }

    # Generate middle values (i.e. not first or last, which are
    # overwritten later)
    if not smooth:
        if dig:
            res = (s >= strike).astype(np.float64)
        else:
            res = np.maximum(s - strike, 0.0)
    else:
        sl = 0.5 * (np.roll(s, 1) + s)
        su = 0.5 * (s + np.roll(s, -1))
        width = su - sl

        if dig:
            res = np.where(
                su <= strike, 0.0, np.where(strike <= sl, 1.0, (su - strike) / width)
            )
        else:
            res = np.where(
                su <= strike,
                0.0,
                np.where(
                    strike <= sl,
                    0.5 * (sl + su) - strike,
                    0.5 * (su - strike) ** 2 / width,
                ),
            )

        res[0] = 1.0 if s[0] >= strike else 0.0 if dig else max(0.0, s[0] - strike)
        res[-1] = 1.0 if s[-1] >= strike else 0.0 if dig else max(0.0, s[-1] - strike)

    if is_put:
        res = 1.0 - res if dig else res - (s - strike)

    return res[None, :]


########################################################################################

# @njit(cache=True, fastmath=True)
def black_scholes_fd(
    spot_price,
    volatility,
    time_to_expiry,
    strike_price,
    risk_free_rate,
    dividend_yield,
    opt_type_value,
    num_steps_per_year=None,
    num_samples=2000,
    num_std=5,
    theta=0.5,
    wind=0,
    digital=False,
    smooth=False,
    update=False,
    return_grid=False,
):

    validate_black_scholes_fd_inputs(
        spot_price,
        volatility,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_yield,
        opt_type_value,
        num_steps_per_year,
        num_samples,
        num_std,
        theta,
        wind,
    )

    if time_to_expiry == 0.0:
        return option_payoff(
            np.array([spot_price]),
            strike_price,
            False,
            digital,
            opt_type_value,
        )[0, 0]

    # Define grid
    std = volatility * (time_to_expiry**0.5)
    xu = num_std * std
    xl = -xu
    d_x = (xu - xl) / max(1, num_samples)
    num_samples = 1 if num_samples <= 0 or xl == xu else num_samples + 1

    # Calculate the drift
    mu = risk_free_rate - dividend_yield

    # Create sample set s
    s = spot_price * np.exp(xl + d_x * np.arange(0, num_samples))

    # Generate the option payoff to be fitted
    payoff = option_payoff(s, strike_price, smooth, digital, opt_type_value)

    # time steps
    if num_steps_per_year is not None:
        num_steps = int(num_steps_per_year * time_to_expiry)
    else:
        num_steps = num_samples // 2

    num_steps = max(num_steps, 1)
    dt = time_to_expiry / max(1, num_steps)

    # Make time series for interest rate, drift, and variance
    r_ = np.zeros(num_samples) + risk_free_rate
    mu_ = mu * s
    var_ = (s * volatility) ** 2

    # Initialise implicit and explicit matricies
    ai = np.empty((num_samples, 3), dtype=np.float64)
    ae = np.empty((num_samples, 3), dtype=np.float64)

    # Store original res as res0
    res = payoff.copy()

    for h in range(num_steps):
        if update or h == 0:
            if theta != 1:
                ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1 - theta, wind)
            if theta != 0:
                ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

        res = fd_roll_backwards(res, theta, ai=ai, ae=ae)

        # Early exercise (American)
        if opt_type_value in {
            OptionTypes.AMERICAN_CALL.value,
            OptionTypes.AMERICAN_PUT.value,
        }:
            idx = res[0] < payoff[0]
            res[0][idx] = payoff[0][idx]

    # By default, keep old behaviour (scalar @ center node)
    center_val = res[0][num_samples // 2]

    # if return_grid:
    #     # Return the full grid for plotting and analysis
    #     return s, res[0], center_val

    return center_val
