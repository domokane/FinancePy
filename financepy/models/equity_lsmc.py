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

@njit(float64(float64, float64, float64, float64, int64, int64, float64, int64,
                 float64, int64, int64), fastmath=True, cache=False)
def equity_lsmc(s0, r, q, v,
                num_steps_per_year,
                num_paths,
                texp,
                option_type_value,
                k,
                use_sobol,
                seed):
    
    if num_paths == 0:
        raise FinError("Num Paths is zero")
    
    np.random.seed(seed)

    num_steps = int(num_steps_per_year * texp)
    num_times = num_steps + 1

    times = np.linspace(0, texp, num_times)
    dt = texp / num_steps
    rootdt = np.sqrt(dt)
    
    vsqrtdt = v * rootdt
    mu = r - q - v*v/2.0
    innov = v * rootdt * np.random.standard_normal((num_times, num_paths))

    if num_paths % 2 == 1:
        num_paths = num_paths + 1

    half_num_paths = int(num_paths/2.0)

    st = np.zeros((num_times, num_paths), 'd')  # stock price matrix
    st[0, :] = s0

    for it in range(1, num_times):
        gp = np.random.standard_normal(half_num_paths)
        g = np.concatenate((gp, -gp))
        st[it, :] = st[it-1, :] * np.exp(mu * dt + v * g * rootdt)

    
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

    value_matrix = np.zeros((exercise_matrix.shape))
    value_matrix[-1] = exercise_matrix[-1]

    # number of basis functions
    poly_degree = 5

    for it in range(num_steps-1, 0, -1):

        df = np.exp(-r * dt)

        regression2 = fit_poly(st[it], 
                               value_matrix[it+1] * df, 
                               poly_degree)

        cont_value = eval_polynomial(regression2, st[it])

        if 1 == 0:
            plt.figure()
            plt.scatter(st[it], value_matrix[it+1])        
            plt.title("it = " + str(it))
            plt.scatter(st[it], cont_value)

        value_matrix[it] = np.where(exercise_matrix[it] > cont_value, 
                                    exercise_matrix[it], 
                                    value_matrix[it+1] * df)

    dt = times[1] - times[0]
    df = np.exp(-r * dt)
    value = df * np.mean(value_matrix[1])
    return value

###############################################################################
