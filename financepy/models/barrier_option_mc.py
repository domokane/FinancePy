# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

from ..utils.error import FinError
from ..utils.global_types import BarrierTypes
from ..models.process_simulator import ProcessSimulator


def value_barrier_option_mc(
    t: float,
    k: float,
    barrier_type: BarrierTypes,
    b: float,
    s: float,
    r: float,
    process_type: int,
    model_params: tuple,
    num_obs_per_year: int = 252,
    num_paths: int = 10000,
    seed: int = 4242,
) -> float:
    """A Monte-Carlo based valuation of the barrier option which simulates
    the evolution of the stock price of at a specified number of annual
    observation times until expiry to examine if the barrier has been
    crossed and the corresponding value of the final payoff, if any. It
    assumes a GBM model for the stock price."""

    if t < 0.0:
        raise FinError("t must be non-negative")

    if k <= 0.0:
        raise FinError("k must be positive")

    if b <= 0.0:
        raise FinError("b must be positive")

    if s <= 0.0:
        raise FinError("s must be positive")

    if num_paths <= 0:
        raise FinError("num_paths must be positive")

    if num_obs_per_year <= 0:
        raise FinError("num_obs_per_year must be positive")

    if not isinstance(barrier_type, BarrierTypes):
        raise FinError(
            "Unknown barrier option type: " + str(barrier_type)
        )

    #######################################################################

    tol = 1.0e-12

    if t == 0.0:

        call_payoff = max(s - k, 0.0)
        put_payoff = max(k - s, 0.0)

        if barrier_type == BarrierTypes.DOWN_AND_OUT_CALL:
            return 0.0 if s <= b + tol else call_payoff

        elif barrier_type == BarrierTypes.DOWN_AND_IN_CALL:
            return call_payoff if s <= b + tol else 0.0

        elif barrier_type == BarrierTypes.UP_AND_OUT_CALL:
            return 0.0 if s >= b - tol else call_payoff

        elif barrier_type == BarrierTypes.UP_AND_IN_CALL:
            return call_payoff if s >= b - tol else 0.0

        elif barrier_type == BarrierTypes.DOWN_AND_OUT_PUT:
            return 0.0 if s <= b + tol else put_payoff

        elif barrier_type == BarrierTypes.DOWN_AND_IN_PUT:
            return put_payoff if s <= b + tol else 0.0

        elif barrier_type == BarrierTypes.UP_AND_OUT_PUT:
            return 0.0 if s >= b - tol else put_payoff

        elif barrier_type == BarrierTypes.UP_AND_IN_PUT:
            return put_payoff if s >= b - tol else 0.0

    # Immediate knock-out
    if barrier_type in (
        BarrierTypes.DOWN_AND_OUT_CALL,
        BarrierTypes.DOWN_AND_OUT_PUT,
    ):
        if s <= b + tol:
            return 0.0

    elif barrier_type in (
        BarrierTypes.UP_AND_OUT_CALL,
        BarrierTypes.UP_AND_OUT_PUT,
    ):
        if s >= b - tol:
            return 0.0

    # Immediate knock-in
    simple_call = (
        barrier_type == BarrierTypes.DOWN_AND_IN_CALL
        and s <= b + tol
    ) or (
        barrier_type == BarrierTypes.UP_AND_IN_CALL
        and s >= b - tol
    )

    simple_put = (
        barrier_type == BarrierTypes.DOWN_AND_IN_PUT
        and s <= b + tol
    ) or (
        barrier_type == BarrierTypes.UP_AND_IN_PUT
        and s >= b - tol
    )

    process = ProcessSimulator()

    # Barrier has already knocked in, so this is now vanilla.
    if simple_call or simple_put:

        s_all = process.get_process(
            process_type,
            t,
            model_params,
            1,
            num_paths,
            seed,
        )

        terminal = s_all[:, -1]

        if simple_call:
            payoff = np.maximum(terminal - k, 0.0)
        else:
            payoff = np.maximum(k - terminal, 0.0)

        return payoff.mean() * np.exp(-r * t)

    # Full barrier simulation.
    num_time_steps = max(
        1,
        int(np.ceil(t * num_obs_per_year)),
    )

    s_all = process.get_process(
        process_type,
        t,
        model_params,
        num_time_steps,
        num_paths,
        seed,
    )

    terminal = s_all[:, -1]

    call_payoff = np.maximum(terminal - k, 0.0)
    put_payoff = np.maximum(k - terminal, 0.0)

    down_crossed = (s_all <= b + tol).any(axis=1)
    up_crossed = (s_all >= b - tol).any(axis=1)

    if barrier_type == BarrierTypes.DOWN_AND_OUT_CALL:
        payoff = call_payoff * ~down_crossed

    elif barrier_type == BarrierTypes.DOWN_AND_IN_CALL:
        payoff = call_payoff * down_crossed

    elif barrier_type == BarrierTypes.UP_AND_OUT_CALL:
        payoff = call_payoff * ~up_crossed

    elif barrier_type == BarrierTypes.UP_AND_IN_CALL:
        payoff = call_payoff * up_crossed

    elif barrier_type == BarrierTypes.DOWN_AND_OUT_PUT:
        payoff = put_payoff * ~down_crossed

    elif barrier_type == BarrierTypes.DOWN_AND_IN_PUT:
        payoff = put_payoff * down_crossed

    elif barrier_type == BarrierTypes.UP_AND_OUT_PUT:
        payoff = put_payoff * ~up_crossed

    elif barrier_type == BarrierTypes.UP_AND_IN_PUT:
        payoff = put_payoff * up_crossed

    else:
        raise FinError(
            "Unknown barrier option type: " + str(barrier_type)
        )

    return payoff.mean() * np.exp(-r * t)

########################################################################################
