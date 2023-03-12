##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..market.curves.discount_curve import DiscountCurve
from ..models.process_simulator import FinProcessSimulator
from ..utils.date import Date
from ..utils.error import FinError
from ..utils.global_types import EquityBarrierTypes
from ..utils.global_vars import gDaysInYear
from ..utils.math import N


def value_one(expiry_date: Date,
              strike_price: float,
              option_type: int,
              barrier_level: float,
              num_observations,  # number of observations per year
              notional: float,
              valuation_date: Date,
              stock_price: (float, np.ndarray),
              rf_rate: float,
              div_rate: float,
              model):
    """ This values a single option. Because of its structure it cannot
    easily be vectorised which is why it has been wrapped. """

    t_exp = (expiry_date - valuation_date) / gDaysInYear

    if t_exp < 0:
        raise FinError("Option expires before value date.")

    t_exp = max(t_exp, 1e-6)

    lnS0k = np.log(stock_price / strike_price)
    sqrtT = np.sqrt(t_exp)

    r = rf_rate
    q = div_rate

    k = strike_price
    s = stock_price
    h = barrier_level

    volatility = model._volatility
    sigma_root_t = volatility * sqrtT
    v2 = volatility * volatility
    mu = r - q
    d1 = (lnS0k + (mu + v2 / 2.0) * t_exp) / sigma_root_t
    d2 = (lnS0k + (mu - v2 / 2.0) * t_exp) / sigma_root_t
    df = np.exp(-r * t_exp)
    dq = np.exp(-q * t_exp)

    c = s * dq * N(d1) - k * df * N(d2)
    p = k * df * N(-d2) - s * dq * N(-d1)

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and s <= h:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and s >= h:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and s >= h:
        return 0.0
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and s <= h:
        return 0.0
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and s <= h:
        return c
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value and s >= h:
        return c
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value and s >= h:
        return p
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and s <= h:
        return p

    num_observations = 1 + t_exp * num_observations

    # Correction by Broadie, Glasserman and Kou, Mathematical Finance, 1997
    # Adjusts the barrier for discrete and not continuous observations
    h_adj = h
    t = t_exp / num_observations

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        h_adj = h * np.exp(0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        h_adj = h * np.exp(-0.5826 * volatility * np.sqrt(t))
    else:
        raise FinError("Unknown barrier option type." +
                       str(option_type))

    h = h_adj

    if abs(volatility) < 1e-5:
        volatility = 1e-5

    l = (mu + v2 / 2.0) / v2
    y = np.log(h * h / (s * k)) / sigma_root_t + l * sigma_root_t
    x1 = np.log(s / h) / sigma_root_t + l * sigma_root_t
    y1 = np.log(h / s) / sigma_root_t + l * sigma_root_t
    hOverS = h / s

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        if h >= k:
            c_do = s * dq * N(x1) - k * df * N(x1 - sigma_root_t) \
                   - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigma_root_t)
            price = c_do
        else:
            c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigma_root_t)
            price = c - c_di
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        if h <= k:
            c_di = s * dq * pow(hOverS, 2.0 * l) * N(y) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * N(y - sigma_root_t)
            price = c_di
        else:
            c_do = s * dq * N(x1) \
                   - k * df * N(x1 - sigma_root_t) \
                   - s * dq * pow(hOverS, 2.0 * l) * N(y1) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * N(y1 - sigma_root_t)
            price = c - c_do
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        if h >= k:
            c_ui = s * dq * N(x1) - k * df * N(x1 - sigma_root_t) \
                   - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   (N(-y + sigma_root_t) - N(-y1 + sigma_root_t))
            price = c_ui
        else:
            price = c
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        if h > k:
            c_ui = s * dq * N(x1) - k * df * N(x1 - sigma_root_t) \
                   - s * dq * pow(hOverS, 2.0 * l) * (N(-y) - N(-y1)) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   (N(-y + sigma_root_t) - N(-y1 + sigma_root_t))
            price = c - c_ui
        else:
            price = 0.0
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        if h > k:
            p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigma_root_t)
            price = p_ui
        else:
            p_uo = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_root_t) \
                   + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   N(-y1 + sigma_root_t)
            price = p - p_uo
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        if h >= k:
            p_ui = -s * dq * pow(hOverS, 2.0 * l) * N(-y) \
                   + k * df * pow(hOverS, 2.0 * l - 2.0) * N(-y + sigma_root_t)
            price = p - p_ui
        else:
            p_uo = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_root_t) \
                   + s * dq * pow(hOverS, 2.0 * l) * N(-y1) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   N(-y1 + sigma_root_t)
            price = p_uo
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        if h >= k:
            price = 0.0
        else:
            p_di = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_root_t) \
                   + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   (N(y - sigma_root_t) - N(y1 - sigma_root_t))
            price = p - p_di
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        if h >= k:
            price = p
        else:
            p_di = -s * dq * N(-x1) \
                   + k * df * N(-x1 + sigma_root_t) \
                   + s * dq * pow(hOverS, 2.0 * l) * (N(y) - N(y1)) \
                   - k * df * pow(hOverS, 2.0 * l - 2.0) * \
                   (N(y - sigma_root_t) - N(y1 - sigma_root_t))
            price = p_di
    else:
        raise FinError("Unknown barrier option type." +
                       str(option_type))

    v = price * notional
    return v


###############################################################################


def value_mc(expiry_date: Date,
             strike_price,
             option_type: int,
             barrier_level,
             notional,
             valuation_date: Date,
             stock_price: float,
             rf_rate: float,
             process_type,
             model_params,
             numAnnObs: int = 252,
             num_paths: int = 10000,
             seed: int = 4242):
    """ A Monte-Carlo based valuation of the barrier option which simulates
    the evolution of the stock price of at a specified number of annual
    observation times until expiry to examine if the barrier has been
    crossed and the corresponding value of the final payoff, if any. It
    assumes a GBM model for the stock price. """

    texp = (expiry_date - valuation_date) / gDaysInYear
    num_time_steps = int(texp * numAnnObs)
    K = strike_price
    B = barrier_level
    option_type = option_type

    process = FinProcessSimulator()

    # TODO - NEED TO DECIDE IF THIS IS PART OF MODEL PARAMS OR NOT ??????????????

    r = rf_rate

    #######################################################################

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value and stock_price <= B:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value and stock_price >= B:
        return 0.0
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value and stock_price <= B:
        return 0.0
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value and stock_price >= B:
        return 0.0

    #######################################################################

    simple_call = False
    simple_put = False

    if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value and stock_price <= B:
        simple_call = True
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value and stock_price >= B:
        simple_call = True
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value and stock_price >= B:
        simple_put = True
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value and stock_price <= B:
        simple_put = True

    if simple_put or simple_call:
        Sall = process.get_process(
            process_type, texp, model_params, 1, num_paths, seed)

    if simple_call:
        c = (np.maximum(Sall[:, -1] - K, 0.0)).mean()
        c = c * np.exp(-r * texp)
        return c

    if simple_put:
        p = (np.maximum(K - Sall[:, -1], 0.0)).mean()
        p = p * np.exp(-r * texp)
        return p

    # Get full set of paths
    Sall = process.get_process(process_type, texp, model_params, num_time_steps,
                               num_paths, seed)

    (num_paths, num_time_steps) = Sall.shape

    if option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value or \
            option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value or \
            option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value or \
            option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:

        barrier_crossed_from_above = [False] * num_paths

        for p in range(0, num_paths):
            barrier_crossed_from_above[p] = np.any(Sall[p] <= B)

    if option_type == EquityBarrierTypes.UP_AND_IN_CALL.value or \
            option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value or \
            option_type == EquityBarrierTypes.UP_AND_IN_PUT.value or \
            option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:

        barrier_crossed_from_below = [False] * num_paths
        for p in range(0, num_paths):
            barrier_crossed_from_below[p] = np.any(Sall[p] >= B)

    payoff = np.zeros(num_paths)
    ones = np.ones(num_paths)

    if option_type == EquityBarrierTypes.DOWN_AND_OUT_CALL.value:
        payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                 (ones - barrier_crossed_from_above)
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_CALL.value:
        payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrier_crossed_from_above
    elif option_type == EquityBarrierTypes.UP_AND_IN_CALL.value:
        payoff = np.maximum(Sall[:, -1] - K, 0.0) * barrier_crossed_from_below
    elif option_type == EquityBarrierTypes.UP_AND_OUT_CALL.value:
        payoff = np.maximum(Sall[:, -1] - K, 0.0) * \
                 (ones - barrier_crossed_from_below)
    elif option_type == EquityBarrierTypes.UP_AND_IN_PUT.value:
        payoff = np.maximum(K - Sall[:, -1], 0.0) * barrier_crossed_from_below
    elif option_type == EquityBarrierTypes.UP_AND_OUT_PUT.value:
        payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                 (ones - barrier_crossed_from_below)
    elif option_type == EquityBarrierTypes.DOWN_AND_OUT_PUT.value:
        payoff = np.maximum(K - Sall[:, -1], 0.0) * \
                 (ones - barrier_crossed_from_above)
    elif option_type == EquityBarrierTypes.DOWN_AND_IN_PUT.value:
        payoff = np.maximum(K - Sall[:, -1], 0.0) * barrier_crossed_from_above
    else:
        raise FinError("Unknown barrier option type." +
                       str(option_type))

    v = payoff.mean() * np.exp(- r * texp)

    return v * notional
