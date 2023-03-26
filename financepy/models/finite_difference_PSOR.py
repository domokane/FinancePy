from copy import deepcopy

import numpy as np

from financepy.utils.math import solve_tridiagonal_matrix, band_matrix_multiplication
from financepy.models.finite_difference import option_payoff
from financepy.utils.global_types import OptionTypes


def calculate_fd_PSOR_matrix(x, risk_free_rate, mu, volatility, dt, theta):
    j = np.arange(len(x))
    alpha = 0.5 * dt * theta * (volatility**2 * j**2 - mu * j)
    beta = 1 - dt * theta * (volatility**2 * j**2 + risk_free_rate)
    kappa = 0.5 * dt * theta * (volatility**2 * j**2 + mu * j)
    return np.array([alpha, beta, kappa])


def black_scholes_fd_PSOR(spot_price, volatility, time_to_expiry,
                          strike_price, risk_free_rate, dividend_yield, option_type,
                          num_steps=None, num_samples=None, s_max=None, theta=0.5, digital=False,
                          smooth=False, delta_bound=1e-5, omega=0.5):
    """
    Solve the Black Scholes equation using the Projected Successive Over Relaxation (PSOR) method.
    This is based on the method outlined here:
    https://www.diva-portal.org/smash/get/diva2:620212/fulltext01.pdf

    This method uses a uniform grid.

    The parameters specific to this model a described below

    Parameters:
        num_steps: Number of time steps used
        num_samples: Number of stock prices the model evalutates
        s_max: Maximum stock price evaluated

        theta: Parameter used for theta solver. Set to 0.5 for Crank Nicolson method

        digital: True when payoff digital
        smooth: True when payoff is smooth

        delta_bound: PSOR with iterate until the square difference over the curve is less than this value
        omega: Weighting parameter for PSOR, balancing previous and current iteration
    """
    if isinstance(option_type, OptionTypes):
        option_type = option_type.value

    # Set default values
    s_max = s_max or max(strike_price, spot_price) * 4
    num_samples = num_samples or s_max * 10
    num_steps = num_steps or int(num_samples // 2)

    # Define grid over stock value and timesteps
    dx = s_max / num_samples
    s = np.arange(s_max, step=dx)
    dt = time_to_expiry / max(1, num_steps)
    timesteps = np.arange(0, num_steps + 1) * dt

    # Define payoff for option as curve in our grid
    payoff = option_payoff(s, strike_price, smooth, digital, option_type)

    # Set drift parameter
    mu = risk_free_rate - dividend_yield

    # Boundary condition is the payoff of the option at expiry.
    # Start here and iterate backwards
    res_k = deepcopy(payoff)[0]

    # Explicit matrix
    Ae = calculate_fd_PSOR_matrix(s, risk_free_rate, mu, volatility, dt, 1 - theta).T
    alpha, beta, kappa = Ae.T

    # Implicit matrix
    Ai = calculate_fd_PSOR_matrix(s, risk_free_rate, mu, volatility, -dt, theta).T
    a, b, c = Ai.T

    # set boundary conditions at s=0 and s=infinity for each timestep
    if option_type in {OptionTypes.EUROPEAN_CALL.value, OptionTypes.AMERICAN_CALL.value}:
        f0 = np.zeros(num_steps + 1)
        fM = s_max - strike_price * np.exp(-risk_free_rate * (time_to_expiry - timesteps))
    elif option_type in {OptionTypes.EUROPEAN_PUT.value, OptionTypes.AMERICAN_PUT.value}:
        f0 = strike_price * np.exp(-risk_free_rate * (time_to_expiry - timesteps))
        fM = np.zeros(num_steps + 1)
    else:
        raise ValueError("Option type not valid for this model")

    for i in range(num_steps - 1, -1, -1):
        # Iterate for until solution is stable
        delta = 1

        # Explicit step
        z_ip1 = band_matrix_multiplication(Ae, 1, 1, res_k)

        # Apply boundary conditions
        z_ip1[0] += alpha[0] * f0[i + 1]
        z_ip1[-1] += kappa[-1] * fM[i + 1]

        # Don't want to change z in place, as we need it later, so deepcopy it
        res_k = deepcopy(z_ip1)

        # Apply boundary conditions
        res_k[0] -= a[0] * f0[i]
        res_k[-1] -= c[-1] * fM[i]

        # Implicit step
        res_k = solve_tridiagonal_matrix(Ai, res_k)

        # PSOR - Iterate until the total change in res is less than delta_bound
        while delta >= delta_bound:
            res_kp1 = (z_ip1 - a * np.roll(res_k, 1) - c * np.roll(res_k, -1)) / b
            res_kp1[0] = res_k[0]
            res_kp1[-1] = res_k[-1]
            res_kp1 = omega * res_k + (1 - omega) * res_kp1

            # Calculate change compared to previous iteration
            delta = np.sum((res_kp1 - res_k) ** 2)

            # roll back i.e. k+1 -> k
            res_k = deepcopy(res_kp1)

        # Early exit for American options
        if option_type in {OptionTypes.AMERICAN_CALL.value, OptionTypes.AMERICAN_PUT.value}:
            idx = res_k < payoff[0]
            res_k[idx] = payoff[0][idx]

    # Interpolate in case spot price doesn't fall exactly on grid value
    return np.interp(spot_price, s, res_k)
