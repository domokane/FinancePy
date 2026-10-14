
import numpy as np
from numba import njit

from ..utils.math import band_matrix_multiplication
from ..utils.global_types import OptionTypes
from .finite_difference import option_payoff, calculate_fd_matrix

########################################################################################


@njit(cache=True, fastmath=True)
def black_scholes_fd_psor(
    spot_price,
    volatility,
    time_to_expiry,
    strike_price,
    risk_free_rate,
    dividend_yield,
    opt_type_value,
    num_time_steps=None,
    num_samples=2000,
    num_std=5,
    theta=0.5,
    wind=0,
    digital=False,
    smooth=False,
    acc=1e-9,
    d_omega=5e-5,
    max_iter=0,
):
    """
    Solve Black-Scholes equation using projected successive over-relaxtion.

    Parameters:
        acc: Keep iterating until this accuracy is achieved
        d_omega: Larger numbers lead to bigger changes in omega with each iteration
        max_iter: Maximum number of iterations in psor step. Set to 0 to
        allow any number of iterations.
    """

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
    num_steps = num_time_steps or num_samples // 4
    dt = time_to_expiry / max(1, num_steps)

    # Make time series for interest rate, drift, and variance
    r_ = np.zeros(num_samples) + risk_free_rate
    mu_ = mu * s
    var_ = (s * volatility) ** 2

    # Result starts as the payoff and we iterate backwards
    # res = deepcopy(payoff)[0]
    res = payoff[0].copy()

    # Explicit
    ae = calculate_fd_matrix(s, r_, mu_, var_, dt, 1.0 - theta, wind)
    # Implicit
    ai = calculate_fd_matrix(s, r_, mu_, var_, -dt, theta, wind)

    # Loop backwards through time_steps, starting with a higher number for
    # faster convergence and lower as estimation improves

    omega = 1.8
    previous_nloops = 0

    total_nloops = 0
    max_nloops = 0

    for _ in range(num_steps - 1, -1, -1):

        res, nloops = psor_roll_backwards(
            res_k=res, acc=acc, ai=ai, ae=ae, omega=omega, max_iter=max_iter
        )

        total_nloops += nloops
        max_nloops = max(max_nloops, nloops)

        # Increase omega if nloops is larger this iteration.
        # Decrease if nloops is less. Otherwise stay the same.
        omega += np.sign(nloops - previous_nloops) * d_omega
        # Omega must be between 1 and 2
        omega = min(max(omega, 1 + d_omega), 2 - d_omega)

        previous_nloops = nloops

        # Early exit for American options
        if opt_type_value in {
            OptionTypes.AMERICAN_CALL.value,
            OptionTypes.AMERICAN_PUT.value,
        }:
            idx = res < payoff[0]
            res[idx] = payoff[0][idx]

    return res[num_samples // 2]


########################################################################################


@njit(fastmath=True, cache=True)
def psor_roll_backwards(ae, ai, res_k, omega, acc=1e-12, max_iter=0):

    # Explicit step
    z_ip1 = band_matrix_multiplication(ae, 1, 1, res_k)

    # psor - Iterate until the total change in res is less than acc
    res_k, nloops = psor(ai, omega, res_k, z_ip1, acc=acc, max_iter=max_iter)
    return res_k, nloops


########################################################################################


@njit(fastmath=True, cache=True)
def psor(ai, omega, res_k, z_ip1, max_iter=0, acc=1e-13):

    nloops = 0
    n = len(res_k)

    a = ai[:, 0]
    b = ai[:, 1]
    c = ai[:, 2]

    while True:

        delta = 0.0

        for j in range(1, n - 1):

            old = res_k[j]

            new = (
                z_ip1[j]
                - a[j] * res_k[j - 1]
                - c[j] * res_k[j + 1]
            ) / b[j]

            new = omega * new + (1.0 - omega) * old

            diff = new - old
            delta += diff * diff

            res_k[j] = new

        nloops += 1

        if delta < acc:
            break

        if max_iter != 0 and nloops >= max_iter:
            raise RuntimeError("Calculation is not converging")

    return res_k, nloops
