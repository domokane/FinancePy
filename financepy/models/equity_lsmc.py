# Copyright (C) 2022 Dominic O'Kane

# https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf

from enum import Enum, auto

import numpy as np

from ..utils.global_types import OptionTypes
from ..utils.polyfit import fit_poly, eval_polynomial
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


# TODO: FIX NUMBA
# @njit(float64(float64, float64, float64, float64, int64, int64, float64, int64,
#              float64, float64, int64, int64, int64), fastmath=True, cache=False)
# @jit

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
    seed,
    ########################################################################################
):

    np.random.seed(seed)

    num_steps = int(num_steps_per_year * time_to_expiry)
    num_times = num_steps + 1

    dt = time_to_expiry / num_times
    times = np.linspace(0, time_to_expiry, num_times)
    rootdt = np.sqrt(dt)

    mu = risk_free_rate - dividend_yield - 0.5 * sigma**2

    if num_paths % 2 == 1:
        num_paths = num_paths + 1

    half_num_paths = int(num_paths / 2.0)

    st = np.zeros((num_times, num_paths), "d")  # stock price matrix
    st[0] = spot_price

    gp = np.random.standard_normal((half_num_paths, num_times))
    g = np.concatenate((gp, -gp))
    for it in range(1, num_times):
        st[it] = st[it - 1] * np.exp(mu * dt + sigma * g[:, it] * rootdt)

    # ensure forward price is recovered exactly
    for it in range(0, num_times):
        fmean = np.mean(st[it])
        fexact = spot_price * np.exp((risk_free_rate - dividend_yield) * times[it])
        st[it] = st[it] * fexact / fmean

    exercise_matrix = np.zeros_like(st)
    for i in range(exercise_matrix.shape[0]):
        exercise_matrix[i] = option_payoff(
            st[i],
            strike_price,
            smooth=False,
            dig=False,
            opt_type=opt_type_value,
        )

    # Set final values for value_matrix and stopping matrix
    value_matrix = np.zeros((exercise_matrix.shape))
    value_matrix[-1] = exercise_matrix[-1]
    stopping = np.zeros_like(value_matrix)
    stopping[-1] = np.where(exercise_matrix[-1] > 0, 1, 0)

    df = np.exp(-risk_free_rate * dt)

    for it in range(num_times - 2, 0, -1):
        if fit_type_value == BoundaryFitTypes.HERMITE_E.value:
            regression2 = np.polynomial.hermite_e.hermefit(
                st[it], value_matrix[it + 1] * df, poly_degree
            )
            cont_value = np.polynomial.hermite_e.hermeval(st[it], regression2)
        elif fit_type_value == BoundaryFitTypes.LAGUERRE.value:
            regression2 = np.polynomial.laguerre.lagfit(
                st[it], value_matrix[it + 1] * df, poly_degree
            )
            cont_value = np.polynomial.laguerre.lagval(st[it], regression2)
        elif fit_type_value == BoundaryFitTypes.HERMITE.value:
            regression2 = np.polynomial.hermite.hermfit(
                st[it], value_matrix[it + 1] * df, poly_degree
            )
            cont_value = np.polynomial.hermite.hermval(st[it], regression2)
        elif fit_type_value == BoundaryFitTypes.LEGENDRE.value:
            regression2 = np.polynomial.legendre.legfit(
                st[it], value_matrix[it + 1] * df, poly_degree
            )
            cont_value = np.polynomial.legendre.legval(st[it], regression2)
        elif fit_type_value == BoundaryFitTypes.CHEBYCHEV.value:
            regression2 = np.polynomial.chebyshev.chebfit(
                st[it], value_matrix[it + 1] * df, poly_degree
            )
            cont_value = np.polynomial.chebyshev.chebval(st[it], regression2)
        elif fit_type_value == BoundaryFitTypes.POLYNOMIAL.value:
            regression2 = fit_poly(st[it], value_matrix[it + 1] * df, poly_degree)
            cont_value = eval_polynomial(regression2, st[it])
        else:
            raise ValueError(f"Unknown _fit_type: {fit_type_value}")
        cont_value[cont_value < 0] = 0

        # Should we exercise at this timestep?
        stopping[it] = np.where(exercise_matrix[it] > cont_value, 1, 0)

        value_matrix[it] = np.where(
            exercise_matrix[it] > cont_value, exercise_matrix[it], cont_value
        )

    # for each path find the earliest stopping time
    values = np.zeros(value_matrix.shape[1])

    for i in range(value_matrix.shape[1]):
        if opt_type_value in {
            OptionTypes.AMERICAN_PUT.value,
            OptionTypes.AMERICAN_CALL.value,
        }:
            # first value in row of stopping matrix that is greater than zero
            s = np.argmax(stopping.T[i])
        else:
            s = -1

        # This is the value of the path discounted to present value
        values[i] = value_matrix[s, i] * np.exp(-risk_free_rate * times[s])

    value = np.mean(values)

    return value
