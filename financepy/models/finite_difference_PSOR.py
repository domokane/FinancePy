from copy import deepcopy

import numpy as np
from numba import njit

from ..utils.math import band_matrix_multiplication
from financepy.utils.global_types import OptionTypes
from .finite_difference import option_payoff, calculate_fd_matrix


def black_scholes_fd_PSOR(spot_price, volatility, time_to_expiry,
                          strike_price, risk_free_rate, dividend_yield, option_type,
                          num_time_steps=None, num_samples=2000, num_std=5, theta=0.5, wind=0, digital=False,
                          smooth=False, acc=1e-13, d_omega=5e-5, max_iter=0):
    """
    Solve Black-Scholes equation using projected successive over-relaxtion.

    Parameters:
        acc: Keep iterating until this accuracy is achieved
        d_omega: Larger numbers lead to bigger changes in omega with each iteration
        max_iter: Maximum number of iterations in PSOR step. Set to 0 to allow any number of iterations.
    """
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

    # Result starts as the payoff and we iterate backwards
    res = deepcopy(payoff)[0]

    # Explicit
    Ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1-theta, wind)
    # Implicit
    Ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

    # Loop backwards through timesteps
    omega = 1.8  # Start with a higher number for faster convergence and lower as estimation improves
    previous_nloops = 0
    for _ in range(num_steps - 1, -1, -1):
        res, nloops = PSOR_roll_backwards(res_k=res, acc=acc, Ai=Ai, Ae=Ae, omega=omega, max_iter=max_iter)

        # Increase omega if nloops is larger this iteration. Decrease if nloops is less. Otherwise stay the same.
        omega += np.sign(nloops - previous_nloops) * d_omega
        omega = min(max(omega, 1 + d_omega), 2 - d_omega)  # Omega must be between 1 and 2

        previous_nloops = nloops

        # Early exit for American options
        if option_type in {OptionTypes.AMERICAN_CALL.value, OptionTypes.AMERICAN_PUT.value}:
            idx = res < payoff[0]
            res[idx] = payoff[0][idx]

    return res[num_samples // 2]


@njit(fastmath=True, cache=True)
def PSOR_roll_backwards(Ae, Ai, res_k, omega, acc=1e-13, max_iter=0):
    # Explicit step
    z_ip1 = band_matrix_multiplication(Ae, 1, 1, res_k)

    # PSOR - Iterate until the total change in res is less than acc
    res_k, nloops = PSOR(Ai, omega, res_k, z_ip1, acc=acc, max_iter=max_iter)
    return res_k, nloops


@njit(fastmath=True, cache=True)
def PSOR(Ai, omega, initial_value, z_ip1, max_iter=0, acc=1e-13):
    """
    Projected Successive Over Relaxation -

    Parameters:
        Ai (np.array): Implicit finite difference matrix
        omega (float): Number between 1 and 2 weighted average of previous and currect iteration
        initial_value (np.array): Vector produced by previous iteration
        z_ip1 (matrix): Matrix for updating to next timestep
        max_iter (int): Maximum number of iterations before raising an error (max_iter=0 means no maximum)
        acc (float): Accuracy. The maximum acceptable square difference between two iterations

    Returns:
        res_k (np.array): Vector for next time step
        nloops (int): The number of iterations required to achieve required accuracy
    """
    nloops = 0
    res_k = initial_value
    res_kp1 = initial_value.copy()
    delta = 1
    a, b, c = Ai.T

    # Iterate for until solution is stable
    while delta >= acc:
        for j in range(1, len(res_k)-1):
            res_kp1[j] = (z_ip1[j] - a[j] * res_kp1[j-1] - c[j] * res_k[j+1]) / b[j]

            # Over relaxation
            res_kp1[j] = omega * res_kp1[j] + (1 - omega) * res_k[j]

        # Calculate change compared to previous iteration
        delta = np.sum((res_kp1 - res_k) ** 2)
        if max_iter != 0 and nloops > max_iter:
            raise RuntimeError("Calculation is not converging")

        # roll back i.e. k+1 -> k
        res_k = res_kp1.copy()

        nloops += 1

    return res_k, nloops
