##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO Fix this

import numpy as np
from numba import njit
from ..finutils.FinGlobalTypes import FinOptionTypes
from ..finutils.FinError import FinError

###############################################################################


@njit(cache=True, fastmath=True)
def _valueMC_NUMBA(s, t, K,  optionType, r, q, v, numPaths, seed):

    # Start pricing here
    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    vsqrtt = v * np.sqrt(t)
    payoff = 0.0

    g = np.random.standard_normal(numPaths)

    if optionType == FinOptionTypes.EUROPEAN_CALL.value:

        for i in range(0, numPaths):
            s_1 = s * np.exp((mu - v2 / 2.0) * t + g[i] * vsqrtt)
            s_2 = s * np.exp((mu - v2 / 2.0) * t - g[i] * vsqrtt)    
            payoff += max(s_1 - K, 0.0)
            payoff += max(s_2 - K, 0.0)

    elif optionType == FinOptionTypes.EUROPEAN_PUT.value:

        for i in range(0, numPaths):
            s_1 = s * np.exp((mu - v2 / 2.0) * t + g[i] * vsqrtt)
            s_2 = s * np.exp((mu - v2 / 2.0) * t - g[i] * vsqrtt)
            payoff += max(K - s_1, 0.0)
            payoff += max(K - s_2, 0.0)

    else:
        raise FinError("Unknown option type.")

    v = payoff * np.exp(-r * t) / numPaths / 2.0
    return v

###############################################################################

def _valueMC_NONUMBA(s, t, K,  optionType, r, q, v, numPaths, seed):

    # Start pricing here
    np.random.seed(seed)
    mu = r - q
    v2 = v**2
    vsqrtt = v * np.sqrt(t)
    payoff = 0.0

    g = np.random.standard_normal(numPaths)

    if optionType == FinOptionTypes.EUROPEAN_CALL.value:

        for i in range(0, numPaths):
            s_1 = s * np.exp((mu - v2 / 2.0) * t + g[i] * vsqrtt)
            s_2 = s * np.exp((mu - v2 / 2.0) * t - g[i] * vsqrtt)    
            payoff += max(s_1 - K, 0.0)
            payoff += max(s_2 - K, 0.0)

    elif optionType == FinOptionTypes.EUROPEAN_PUT.value:

        for i in range(0, numPaths):
            s_1 = s * np.exp((mu - v2 / 2.0) * t + g[i] * vsqrtt)
            s_2 = s * np.exp((mu - v2 / 2.0) * t - g[i] * vsqrtt)
            payoff += max(K - s_1, 0.0)
            payoff += max(K - s_2, 0.0)

    else:
        raise FinError("Unknown option type.")

    v = payoff * np.exp(-r * t) / numPaths / 2.0
    return v

###############################################################################
