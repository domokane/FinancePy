from copy import deepcopy

import numpy as np
from numba import njit

from ..utils.math import band_matrix_multiplication
from financepy.utils.global_types import OptionTypes
from .finite_difference import option_payoff, calculate_fd_matrix


def black_scholes_fd_PSOR(spot_price, volatility, time_to_expiry,
                          strike_price, risk_free_rate, dividend_yield, option_type,
                          num_time_steps=None, num_samples=5000, num_std=5, theta=0.5, wind=0, digital=False,
                          smooth=False):
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
    num_steps = num_time_steps or num_samples // 4
    dt = time_to_expiry / max(1, num_steps)

    # Make time series for interest rate, drift, and variance
    r_ = np.zeros(num_samples) + risk_free_rate
    mu_ = mu * s
    var_ = (s * volatility) ** 2

    # Store original res as res0
    res = deepcopy(payoff)[0]

    # Explicit
    Ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1-theta, wind)
    # Implicit
    Ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

    res = PSOR_roll_backwards(res_k=res, payoff=payoff, option_type=option_type, delta_bound=1e-5, Ai=Ai, Ae=Ae,
                              num_steps=num_steps)

    return res[num_samples // 2]


@njit(fastmath=True, cache=True)
def PSOR_roll_backwards(Ae, Ai, delta_bound, option_type, payoff, res_k, num_steps):
    # Loop backwards through timesteps
    for i in range(num_steps - 1, -1, -1):
        # Explicit step
        z_ip1 = band_matrix_multiplication(Ae, 1, 1, res_k)

        # PSOR - Iterate until the total change in res is less than delta_bound
        omega = 1.9  # TODO Dynamically update omega
        res_k = PSOR(Ai, delta_bound, omega, res_k, z_ip1)

        # Early exit for American options
        if option_type in {OptionTypes.AMERICAN_CALL.value, OptionTypes.AMERICAN_PUT.value}:
            idx = res_k < payoff[0]
            res_k[idx] = payoff[0][idx]
    return res_k


@njit(fastmath=True, cache=True)
def PSOR(Ai, acc, omega, initial_value, z_ip1, max_iter=1000):
    nloops = 0
    res_k = initial_value
    res_kp1 = initial_value.copy()
    delta = 1
    a, b, c = Ai.T

    # Iterate for until solution is stable
    while delta >= acc:
        for j in range(1, len(res_k)-1):
            res_kp1[j] = (z_ip1[j] - a[j] * res_kp1[j-1] - c[j] * res_k[j+1]) / b[j]

            res_kp1[j] = omega * res_kp1[j] + (1 - omega) * res_k[j]

        # Calculate change compared to previous iteration
        delta = np.sum((res_kp1 - res_k) ** 2)
        if nloops > max_iter:
            raise RuntimeError("Calculation is not converging")

        # roll back i.e. k+1 -> k
        res_k = res_kp1.copy()

        nloops +=1

    return res_k
