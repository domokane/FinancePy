# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from ..utils.error import FinError
from ..utils.global_types import EquityBarrierTypes

from ..models.process_simulator import FinProcessSimulator
from ..models.process_simulator import ProcessTypes


def value_equity_barrier_option_mc(
    t: float,
    k: float,
    opt_type: EquityBarrierTypes,
    b: float,
    notional: float,
    s: float,
    r: float,
    process_type: int,
    model_params: dict,
    num_obs_per_year: int = 252,
    num_paths: int = 10000,
    seed: int = 4242,
) -> float:
    """A Monte-Carlo based valuation of the barrier option which simulates
    the evolution of the stock price of at a specified number of annual
    observation times until expiry to examine if the barrier has been
    crossed and the corresponding value of the final payoff, if any. It
    assumes a GBM model for the stock price."""

    tol = 1e-12  # barrier tolerance

    t = max(t, 1e-6)
    num_time_steps = int(t * num_obs_per_year)

    process_type = ProcessTypes.GBM_PROCESS

    process = FinProcessSimulator()

    #######################################################################

    if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL and s <= b + tol:
        return 0.0
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL and s >= b - tol:
        return 0.0
    elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT and s <= b + tol:
        return 0.0
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT and s >= b - tol:
        return 0.0

    #######################################################################

    simple_call = False
    simple_put = False

    if opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL and s <= b + tol:
        simple_call = True
    elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL and s >= b - tol:
        simple_call = True
    elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT and s >= b + tol:
        simple_put = True
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT and s <= b - tol:
        simple_put = True

    s_all = None

    # For a simple call or put we only need the terminal stock price distribution
    if simple_put or simple_call:
        s_all = process.get_process(process_type, t, model_params, 1, num_paths, seed)

    if simple_call:
        c = (np.maximum(s_all[:, -1] - k, 0.0)).mean()
        c = c * np.exp(-r * t)
        return c

    if simple_put:
        p = (np.maximum(k - s_all[:, -1], 0.0)).mean()
        p = p * np.exp(-r * t)
        return p

    # Get full set of paths
    s_all = process.get_process(
        process_type, t, model_params, num_time_steps, num_paths, seed
    )

    (num_paths, num_time_steps) = s_all.shape

    if (
        opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL
        or opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL
        or opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT
        or opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT
    ):

        barrier_crossed_from_above = [False] * num_paths

        for p in range(0, num_paths):
            barrier_crossed_from_above[p] = np.any(s_all[p] <= b)

    if (
        opt_type == EquityBarrierTypes.UP_AND_IN_CALL
        or opt_type == EquityBarrierTypes.UP_AND_OUT_CALL
        or opt_type == EquityBarrierTypes.UP_AND_IN_PUT
        or opt_type == EquityBarrierTypes.UP_AND_OUT_PUT
    ):

        barrier_crossed_from_below = [False] * num_paths
        for p in range(0, num_paths):
            barrier_crossed_from_below[p] = np.any(s_all[p] >= b)

    payoff = np.zeros(num_paths)
    ones = np.ones(num_paths)

    if opt_type == EquityBarrierTypes.DOWN_AND_OUT_CALL:
        crossed = (s_all <= b + tol).any(axis=1)
        payoff = np.maximum(s_all[:, -1] - k, 0.0) * (ones - crossed)
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_CALL:
        crossed = (s_all <= b + tol).any(axis=1)
        payoff = np.maximum(s_all[:, -1] - k, 0.0) * crossed
    elif opt_type == EquityBarrierTypes.UP_AND_IN_CALL:
        crossed = (s_all >= b - tol).any(axis=1)
        payoff = np.maximum(s_all[:, -1] - k, 0.0) * crossed
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_CALL:
        crossed = (s_all >= b - tol).any(axis=1)
        payoff = np.maximum(s_all[:, -1] - k, 0.0) * (ones - crossed)
    elif opt_type == EquityBarrierTypes.UP_AND_IN_PUT:
        crossed = (s_all >= b - tol).any(axis=1)
        payoff = np.maximum(k - s_all[:, -1], 0.0) * crossed
    elif opt_type == EquityBarrierTypes.UP_AND_OUT_PUT:
        crossed = (s_all >= b - tol).any(axis=1)
        payoff = np.maximum(k - s_all[:, -1], 0.0) * (ones - crossed)
    elif opt_type == EquityBarrierTypes.DOWN_AND_OUT_PUT:
        crossed = (s_all <= b + tol).any(axis=1)
        payoff = np.maximum(k - s_all[:, -1], 0.0) * (ones - crossed)
    elif opt_type == EquityBarrierTypes.DOWN_AND_IN_PUT:
        crossed = (s_all <= b + tol).any(axis=1)
        payoff = np.maximum(k - s_all[:, -1], 0.0) * crossed
    else:
        raise FinError("Unknown barrier option type." + str(opt_type))

    v = payoff.mean() * np.exp(-r * t)

    return v * notional


########################################################################################
