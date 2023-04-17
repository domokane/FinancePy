###############################################################################
# Copyright (C) 2022 Dominic O'Kane
###############################################################################


# https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf

import sys

import numpy as np
from numba import jit, njit, float64, int64

from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from ..utils.polyfit import fit_poly, eval_polynomial
from ..models.sobol import get_gaussian_sobol

import matplotlib.pyplot as plt

# This is a first implementation of American Monte Carlo using the method of 
# Longstaff and Schwartz. Work is needed to add laguerre Polynomials and 
# other interpolation methods.

#@njit(float64(float64, float64, float64, float64, int64, int64, float64, int64,
#                 float64, int64, int64, int64), fastmath=True, cache=False)
def equity_lsmc(s0, r, q, v,
                num_steps_per_year,
                num_paths,
                texp,
                option_type_value,
                k,
                poly_degree,
                use_sobol,
                seed):
    
    if num_paths == 0:
        raise FinError("Num Paths is zero")
    
    np.random.seed(seed)

    num_steps = int(num_steps_per_year * texp)
    num_times = num_steps + 1

    dt = texp / num_times
    times = np.linspace(0, texp, num_times)
    rootdt = np.sqrt(dt)
    
    mu = r - q - 0.5 * v**2

    if num_paths % 2 == 1:
        num_paths = num_paths + 1

    half_num_paths = int(num_paths/2.0)

    st = np.zeros((num_times, num_paths), 'd')  # stock price matrix
    st[0] = s0

    for it in range(1, num_times):
        gp = np.random.standard_normal(half_num_paths)
        g = np.concatenate((gp, -gp))
        st[it] = st[it-1] * np.exp(mu * dt + v * g * rootdt)

    
    # ensure forward price is recovered exactly
    for it in range(0, num_times):
        t = times[it]
        fmean = np.mean(st[it])
        fexact = s0 * np.exp((r-q)*t)
        st[it] = st[it] * fexact / fmean
        
    if option_type_value == OptionTypes.AMERICAN_CALL.value:
        exercise_matrix = np.maximum(st - k, 0.0)
    elif option_type_value == OptionTypes.AMERICAN_PUT.value:
        exercise_matrix = np.maximum(k - st, 0.0)
    else:
        raise FinError("Unknown option type.")

    # Set final values for value_matrix and stopping matrix
    value_matrix = np.zeros((exercise_matrix.shape))
    value_matrix[-1] = exercise_matrix[-1]
    stopping = np.zeros_like(value_matrix)
    stopping[-1] = np.where(exercise_matrix[-1] > 0, 1, 0)

    df = np.exp(-r * dt)
    for it in range(num_times-2, 0, -1):
        regression2 = np.polynomial.laguerre.lagfit(st[it], value_matrix[it + 1] * df, poly_degree)
        cont_value = np.polynomial.laguerre.lagval(st[it], regression2)
        #vander = np.vander(st[it], poly_degree+1)
        #regression2 = np.linalg.lstsq(vander, value_matrix[it+1] * df)[0]
        #cont_value = eval_polynomial(regression2, st[it])
        cont_value[cont_value < 0] = 0

        # Should we exercise at this timestep?
        stopping[it] = np.where(exercise_matrix[it] > cont_value, 1, 0)

        value_matrix[it] = np.where(exercise_matrix[it] > cont_value,
                                    exercise_matrix[it], 
                                    cont_value)

    # for each path find the earliest stopping time
    values = np.zeros(value_matrix.shape[1])
    for i in range(value_matrix.shape[1]):
        # first value in row of stopping matrix that is greater than zero
        s = np.argmax(stopping.T[i])

        # This is the value of the path discounted to present value
        values[i] = value_matrix[s, i] * np.exp(-r * times[s])
    value = np.mean(values)

    return value

###############################################################################
