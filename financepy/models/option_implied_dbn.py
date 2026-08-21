# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

# TODO Fix this

import numpy as np

from ..utils.global_types import OptionTypes
from ..utils.error import FinError

from .black_scholes_analytic import european_value
from typing import Sequence

###############################################################################


def option_implied_dbn(
    s: float,
    t: float,
    r: float,
    q: float,
    strikes: Sequence[float],
    sigmas: Sequence[float],
) -> np.ndarray:
    """Calculate the option smile/skew-implied risk-neutral distribution.

    Uses the Breeden-Litzenberger result

        f(K) = exp(rT) * d²C/dK²

    and returns

        f(K) * dK

    at each strike, rather than the raw probability density f(K).

    Parameters
    ----------
    s : float
        Current underlying spot price.
    t : float
        Time to expiry in years.
    r : float
        Continuously compounded risk-free rate.
    q : float
        Continuously compounded dividend yield.
    strikes : Sequence[float]
        Increasing, equally spaced strike grid.
    sigmas : Sequence[float]
        Black-Scholes implied volatility corresponding to each strike.

    Returns
    -------
    np.ndarray
        Probability mass f(K) * dK associated with each strike-grid point.
        The first and last values are zero because a central finite
        difference cannot be calculated there.
    """

    strikes = np.asarray(strikes, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)

    if strikes.ndim != 1 or sigmas.ndim != 1:
        raise FinError("Strikes and sigmas must be one-dimensional.")

    if len(strikes) != len(sigmas):
        raise FinError("Strike and Sigma vectors do not have same length.")

    num_steps = len(strikes)

    if num_steps < 3:
        raise FinError("At least three strikes are required.")

    if s <= 0.0:
        raise FinError("Spot price must be positive.")

    if t <= 0.0:
        raise FinError("Time to expiry must be positive.")

    if np.any(strikes <= 0.0):
        raise FinError("Strikes must be positive.")

    if np.any(sigmas <= 0.0):
        raise FinError("Volatilities must be positive.")

    strike_diffs = np.diff(strikes)

    if np.any(strike_diffs <= 0.0):
        raise FinError("Strikes must be strictly increasing.")

    # This finite-difference implementation assumes an equally spaced grid.
    dk = strike_diffs[0]

    if not np.allclose(strike_diffs, dk):
        raise FinError("Strikes must be equally spaced.")

    # Calculate European call values at every strike using the
    # strike-dependent implied volatility.
    values = np.empty(num_steps)

    for ik in range(num_steps):
        values[ik] = european_value(
            s,
            t,
            strikes[ik],
            r,
            q,
            sigmas[ik],
            OptionTypes.EUROPEAN_CALL.value,
        )

    # Breeden-Litzenberger:
    #
    #     f(K) = exp(rT) * d²C/dK²
    #
    # Central finite difference:
    #
    #              C(K+dK) - 2C(K) + C(K-dK)
    #     d²C/dK² ≈ --------------------------
    #                         dK²
    #
    # This function returns f(K) * dK, so:
    #
    #     f(K) dK ≈ exp(rT)
    #                * [C(K+dK) - 2C(K) + C(K-dK)] / dK
    #

    density_dk = np.zeros(num_steps)

    inflator = np.exp(r * t)

    second_diffs = values[2:] - 2.0 * values[1:-1] + values[:-2]

    density_dk[1:-1] = inflator * second_diffs / dk

    return density_dk


########################################################################################