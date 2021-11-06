##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from numba import njit, float64

from ..utils.global_types import OptionTypes
from ..utils.error import FinError

from .black_scholes_analytic import bs_value

###############################################################################
# Analytical Black Scholes model implementation and approximations
###############################################################################


@njit(float64[:](float64, float64, float64, float64, float64[:],
                 float64[:]), cache=True, fastmath=True)
def option_implied_dbn(s, t, r, q, strikes, sigmas):
    """ This function calculates the option smile/skew-implied probability
    density function times the interval width. """

    if len(strikes) != len(sigmas):
        raise FinError("Strike and Sigma vector do not have same length.")

    num_steps = len(strikes)

    sigma = sigmas[0]
    strike = strikes[0]

    sigma = sigmas[1]
    strike = strikes[1]

    inflator = np.exp((r-0) * t)
    dK = strikes[1] - strikes[0]
    values = np.zeros(num_steps)

    for ik in range(0, num_steps):
        strike = strikes[ik]
        sigma = sigmas[ik]
        v = bs_value(s, t, strike, r, q, sigma,
                     OptionTypes.EUROPEAN_CALL.value)
        values[ik] = v

    # Calculate the density rho(K) dK
    densitydk = np.zeros(num_steps)

    for ik in range(1, num_steps-1):
        d2VdK2 = (values[ik+1] - 2.0 * values[ik] + values[ik-1]) / dK

 #       print("%d %12.8f %12.8f %12.8f" %
 #             (ik, strikes[ik], values[ik], d2VdK2))

        densitydk[ik] = d2VdK2 * inflator

    return densitydk

###############################################################################
