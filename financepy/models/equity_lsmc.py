# Copyright (C) 2022 Dominic O'Kane

# https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf

from enum import Enum, auto

import numpy as np

from ..utils.global_types import OptionTypes
from ..models.finite_difference import option_payoff

# This is a first implementation of American Monte Carlo using the method of
# Longstaff and Schwartz. Work is needed to add laguerre Polynomials and
# other interpolation methods.

########################################################################################


class BoundaryFitTypes(Enum):
    """Enum for polynomial fit types."""

    HERMITE_E = auto()
    LAGUERRE = auto()
    HERMITE = auto()
    LEGENDRE = auto()
    CHEBYCHEV = auto()
    POLYNOMIAL = auto()


########################################################################################
# Many functions are not numba compatible so no benefit
# @njit(cache=True, fastmath=True)
########################################################################################

def equity_lsmc(
    spot_price,
    risk_free_rate,
    dividend_yield,
    sigma,
    num_paths,
    num_steps_per_year,
    time_to_expiry,
    opt_type_value,
    strike_price,
    poly_degree,
    fit_type_value,
    use_sobol,
    seed,  # Unused
):

    np.random.seed(seed)

    num_steps = int(num_steps_per_year * time_to_expiry)
    num_steps = max(num_steps, 1)

    num_times = num_steps + 1
    dt = time_to_expiry / num_steps
    times = np.linspace(0.0, time_to_expiry, num_times)
    rootdt = np.sqrt(dt)

    mu = risk_free_rate - dividend_yield - 0.5 * sigma**2

    if num_paths % 2 == 1:
        num_paths += 1

    half_num_paths = num_paths // 2

    st = np.zeros((num_times, num_paths))
    st[0] = spot_price

    gp = np.random.standard_normal((half_num_paths, num_times - 1))
    g = np.concatenate((gp, -gp))

    for it in range(1, num_times):
        st[it] = st[it - 1] * np.exp(mu * dt + sigma * g[:, it - 1] * rootdt)

    # ensure forward price is recovered exactly
    for it in range(0, num_times):
        fmean = np.mean(st[it])
        fexact = spot_price * np.exp((risk_free_rate - dividend_yield) * times[it])
        st[it] = st[it] * fexact / fmean

    exercise_matrix = np.zeros_like(st)
    for it in range(0, num_times):
        exercise_matrix[it] = option_payoff(
            st[it],
            strike_price,
            smooth=False,
            dig=False,
            opt_type_int=opt_type_value,
        )

    # Set final values for value_matrix and stopping matrix
    value_matrix = np.zeros((exercise_matrix.shape))
    value_matrix[-1] = exercise_matrix[-1]

    stopping = np.zeros_like(value_matrix)
    stopping[-1] = np.where(exercise_matrix[-1] > 0, 1, 0)

    df = np.exp(-risk_free_rate * dt)

    for it in range(num_times - 2, 0, -1):

        itm = exercise_matrix[it] > 0.0

        if np.sum(itm) <= poly_degree:
            value_matrix[it] = value_matrix[it + 1] * df
            continue

        else:
            x = st[it][itm]
            y = value_matrix[it + 1][itm] * df

        if fit_type_value == BoundaryFitTypes.HERMITE_E.value:
            regression = np.polynomial.HermiteE.fit(x, y, poly_degree)

        elif fit_type_value == BoundaryFitTypes.LAGUERRE.value:
            regression = np.polynomial.Laguerre.fit(x, y, poly_degree)

        elif fit_type_value == BoundaryFitTypes.HERMITE.value:
            regression = np.polynomial.Hermite.fit(x, y, poly_degree)

        elif fit_type_value == BoundaryFitTypes.LEGENDRE.value:
            regression = np.polynomial.Legendre.fit(x, y, poly_degree)

        elif fit_type_value == BoundaryFitTypes.CHEBYCHEV.value:
            regression = np.polynomial.Chebyshev.fit(x, y, poly_degree)

        elif fit_type_value == BoundaryFitTypes.POLYNOMIAL.value:
            regression = np.polynomial.Polynomial.fit(x, y, poly_degree)

        else:
            raise ValueError(f"Unknown fit type: {fit_type_value}")

        cont_value = np.maximum(regression(st[it]), 0.0)

        exercise_now = itm & (exercise_matrix[it] > cont_value)

        stopping[it] = np.where(exercise_now, 1.0, 0.0)

        value_matrix[it] = np.where(
            exercise_now, exercise_matrix[it], value_matrix[it + 1] * df
        )

    # for each path find the earliest stopping time
    values = np.zeros(value_matrix.shape[1])

    is_american = opt_type_value in {
        OptionTypes.AMERICAN_CALL.value,
        OptionTypes.AMERICAN_PUT.value,
    }

    for i in range(0, num_paths):

        if is_american:
            exercise_times = np.where(stopping[:, i] > 0.0)[0]

            if len(exercise_times) > 0:
                stop_time = exercise_times[0]
            else:
                stop_time = num_times - 1
        else:
            stop_time = num_times - 1

        # This is the value of the path discounted to present value
        values[i] = value_matrix[stop_time, i] * np.exp(
            -risk_free_rate * times[stop_time]
        )

    value = np.mean(values)

    if is_american:
        value = max(value, exercise_matrix[0, 0])

    return value
