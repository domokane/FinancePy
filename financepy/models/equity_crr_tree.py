# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


from typing import Dict, Any
import numpy as np
from numba import njit

from ..utils.global_types import OptionTypes
from ..utils.error import FinError


########################################################################################


@njit(fastmath=True, cache=True)
def crr_tree_val(
    stock_price,
    interest_rate,
    dividend_rate,
    volatility,
    num_steps_per_year,
    time_to_expiry,
    opt_type,
    strike_price,
    is_even,
) -> np.ndarray:
    """
    Value an American option using a Binomial Tree.

    Args:
        stock_price (float): Initial stock price.
        interest_rate (float): Continuously compounded risk-free rate.
        dividend_rate (float): Continuously compounded dividend yield.
        volatility (float): Black-Scholes volatility.
        num_steps_per_year (int): Number of tree steps per year.
        time_to_expiry (float): Time to expiry in years.
        opt_type (int): Option type (OptionTypes enum value).
        strike_price (float): Option strike price.
        is_even (int): 1 for even, 0 for odd number of steps.
    Returns:
        np.ndarray: Array of [price, delta, gamma, theta].
    """

    num_steps = int(num_steps_per_year * time_to_expiry)
    num_steps = max(num_steps, 30)

    # if the number of steps is even but we want odd then make it odd
    if num_steps % 2 == 0 and is_even == 0:
        num_steps += 1
    elif num_steps % 2 == 1 and is_even == 1:
        num_steps += 1

    is_call = (
        opt_type == OptionTypes.EUROPEAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_CALL.value
    )

    is_american = (
        opt_type == OptionTypes.AMERICAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_PUT.value
    )

    #    print(num_steps)
    # this is the size of the step
    dt = time_to_expiry / num_steps
    r = interest_rate
    q = dividend_rate

    # the number of nodes on the tree
    num_nodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    stock_values = np.zeros(num_nodes)
    stock_values[0] = stock_price

    option_values = np.zeros(num_nodes)
    u = np.exp(volatility * np.sqrt(dt))
    d = 1.0 / u
    s_low = stock_price

    a = np.exp((r - q) * dt)
    p = (a - d) / (u - d)
    df = np.exp(-r * dt)

    for i_time in range(1, num_steps + 1):
        s_low *= d
        s = s_low
        for i_node in range(0, i_time + 1):
            index = i_time * (i_time + 1) // 2
            stock_values[index + i_node] = s
            s = s * (u * u)

    # work backwards by first setting values at expiry date
    index = int(0.5 * num_steps * (num_steps + 1))

    for i_node in range(0, num_steps + 1):

        s = stock_values[index + i_node]

        if is_call:
            option_values[index + i_node] = np.maximum(s - strike_price, 0.0)
        else:
            option_values[index + i_node] = np.maximum(strike_price - s, 0.0)

    # begin backward steps from expiry to value date
    for i_time in range(num_steps - 1, -1, -1):

        index = i_time * (i_time + 1) // 2

        for i_node in range(0, i_time + 1):

            s = stock_values[index + i_node]

            exercise_value = 0.0

            if not is_american:
                exercise_value = 0.0
            elif is_call:
                exercise_value = np.maximum(s - strike_price, 0.0)
            else:
                exercise_value = np.maximum(strike_price - s, 0.0)

            next_index = (i_time + 1) * (i_time + 2) // 2

            next_node_dn = next_index + i_node
            next_node_up = next_index + i_node + 1

            v_up = option_values[next_node_up]
            v_dn = option_values[next_node_dn]
            future_expected_value = p * v_up
            future_expected_value += (1.0 - p) * v_dn
            hold_value = df * future_expected_value

            if is_american is False:
                option_values[index + i_node] = hold_value
            else:
                option_values[index + i_node] = np.maximum(exercise_value, hold_value)

    # We calculate all of the important Greeks in one go
    price = option_values[0]
    delta = (option_values[2] - option_values[1]) / (stock_values[2] - stock_values[1])
    delta_up = (option_values[5] - option_values[4]) / (
        stock_values[5] - stock_values[4]
    )
    delta_dn = (option_values[4] - option_values[3]) / (
        stock_values[4] - stock_values[3]
    )
    gamma = (delta_up - delta_dn) / (stock_values[2] - stock_values[1])
    theta = (option_values[4] - option_values[0]) / (2.0 * dt)
    results = np.array([price, delta, gamma, theta])
    return results

#########################################################################################


@njit(cache=True, fastmath=True)
def crr_tree_val_fast(
    stock_price,
    interest_rate,
    dividend_rate,
    volatility,
    total_num_steps,  # We do not do steps per year here
    time_to_expiry,
    opt_type,
    strike_price,
):

    num_steps = max(total_num_steps, 10)

    is_call = (
        opt_type == OptionTypes.EUROPEAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_CALL.value
    )
    is_american = (
        opt_type == OptionTypes.AMERICAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_PUT.value
    )

    dt = time_to_expiry / num_steps
    u = np.exp(volatility * np.sqrt(dt))
    d = 1.0 / u

    a = np.exp((interest_rate - dividend_rate) * dt)
    p = (a - d) / (u - d)

    if p < 0.0 or p > 1.0:
        raise FinError("Invalid CRR probability")

    df = np.exp(-interest_rate * dt)

    # Terminal stock values
    values = np.empty(num_steps + 1)

    s_low = stock_price * (d**num_steps)

    s = s_low
    uu = u * u

    for j in range(num_steps + 1):
        if is_call:
            values[j] = max(s - strike_price, 0.0)
        else:
            values[j] = max(strike_price - s, 0.0)
        s *= uu

    # Need these nodes later for Greeks.
    v10 = 0.0
    v11 = 0.0
    v20 = 0.0
    v21 = 0.0
    v22 = 0.0

    # Roll backwards
    for i in range(num_steps - 1, -1, -1):

        # Lowest stock at level i
        s_low *= u
        s = s_low

        for j in range(i + 1):
            hold = df * (
                p * values[j + 1]
                + (1.0 - p) * values[j]
            )

            if is_american:
                if is_call:
                    exercise = max(s - strike_price, 0.0)
                else:
                    exercise = max(strike_price - s, 0.0)

                values[j] = max(exercise, hold)
            else:
                values[j] = hold

            s *= uu

        # Preserve levels required for Greeks
        if i == 2:
            v20 = values[0]
            v21 = values[1]
            v22 = values[2]
        elif i == 1:
            v10 = values[0]
            v11 = values[1]

    price = values[0]

    s_dn = stock_price * d
    s_up = stock_price * u

    delta = (v11 - v10) / (s_up - s_dn)

    s_dd = stock_price * d * d
    s_ud = stock_price
    s_uu = stock_price * u * u

    delta_dn = (v21 - v20) / (s_ud - s_dd)
    delta_up = (v22 - v21) / (s_uu - s_ud)

    gamma = (delta_up - delta_dn) / (s_up - s_dn)
    theta = (v21 - price) / (2.0 * dt)

    return np.array([price, delta, gamma, theta])


########################################################################################


@njit(cache=True, fastmath=True)
def crr_tree_val_avg(
    stock_price: float,
    interest_rate: float,
    dividend_rate: float,
    volatility: float,
    num_steps_per_year: int,
    time_to_expiry: float,
    opt_type: int,
    strike_price: float,
    do_averaging: bool,
) -> Dict[str, Any]:
    """Calculate CRR tree values, optionally averaging adjacent
    even and odd trees."""

    num_steps = int(num_steps_per_year * time_to_expiry)

    value1 = crr_tree_val_fast(
        stock_price,
        interest_rate,
        dividend_rate,
        volatility,
        num_steps,
        time_to_expiry,
        opt_type,
        strike_price,
    )

    if do_averaging is False:
        v = value1

    else:

        value2 = crr_tree_val_fast(
            stock_price,
            interest_rate,
            dividend_rate,
            volatility,
            num_steps+1,
            time_to_expiry,
            opt_type,
            strike_price,
        )

        v = (value1 + value2) / 2.0

    res = {
        "value": v[0],
        "delta": v[1],
        "gamma": v[2],
        "theta": v[3],
    }

    return res

#########################################################################################


def crr_richardson(
    stock_price,
    interest_rate,
    dividend_rate,
    volatility,
    num_steps_per_year,
    time_to_expiry,
    opt_type,
    strike_price,
):

    num_steps = int(num_steps_per_year * time_to_expiry)

    # N and N+1
    v_n = crr_tree_val_fast(
        stock_price,
        interest_rate,
        dividend_rate,
        volatility,
        num_steps,
        time_to_expiry,
        opt_type,
        strike_price,
    )

    v_n1 = crr_tree_val_fast(
        stock_price,
        interest_rate,
        dividend_rate,
        volatility,
        num_steps + 1,
        time_to_expiry,
        opt_type,
        strike_price,
    )

    v_avg_n = 0.5 * (v_n + v_n1)

    # 2N and 2N+1
    v_2n = crr_tree_val_fast(
        stock_price,
        interest_rate,
        dividend_rate,
        volatility,
        2 * num_steps,
        time_to_expiry,
        opt_type,
        strike_price,
    )

    v_2n1 = crr_tree_val_fast(
        stock_price,
        interest_rate,
        dividend_rate,
        volatility,
        2 * num_steps + 1,
        time_to_expiry,
        opt_type,
        strike_price,
    )

    v_avg_2n = 0.5 * (v_2n + v_2n1)

    # Richardson extrapolation
    v = 2.0 * v_avg_2n - v_avg_n

    return {
        "value": v[0],
        "delta": v[1],
        "gamma": v[2],
        "theta": v[3],
    }


########################################################################################

@njit(cache=True, fastmath=True)
def _pp2_inversion(z, n):
    """Peizer-Pratt inversion used by the Leisen-Reimer tree."""

    denom = n + 1.0 / 3.0 + 0.1 / (n + 1.0)

    x = (
        (z / denom) ** 2
        * (n + 1.0 / 6.0)
    )

    if z >= 0.0:
        return 0.5 + 0.5 * np.sqrt(1.0 - np.exp(-x))
    else:
        return 0.5 - 0.5 * np.sqrt(1.0 - np.exp(-x))

########################################################################################


@njit(cache=True, fastmath=True)
def leisen_reimer_tree(
    stock_price,
    interest_rate,
    dividend_rate,
    volatility,
    num_steps,
    time_to_expiry,
    opt_type,
    strike_price,
):

    # LR is normally constructed with an odd number of steps.
    num_steps = max(num_steps, 11)

    if num_steps % 2 == 0:
        num_steps += 1

    is_call = (
        opt_type == OptionTypes.EUROPEAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_CALL.value
    )

    is_american = (
        opt_type == OptionTypes.AMERICAN_CALL.value
        or opt_type == OptionTypes.AMERICAN_PUT.value
    )

    dt = time_to_expiry / num_steps

    r = interest_rate
    q = dividend_rate
    sigma = volatility
    T = time_to_expiry

    sqrt_t = np.sqrt(T)

    d1 = (
        np.log(stock_price / strike_price)
        + (r - q + 0.5 * sigma * sigma) * T
    ) / (sigma * sqrt_t)

    d2 = d1 - sigma * sqrt_t

    # Leisen-Reimer probabilities.
    p = _pp2_inversion(d2, num_steps)
    p_prime = _pp2_inversion(d1, num_steps)

    # Expected growth over one step.
    a = np.exp((r - q) * dt)

    # LR up/down factors.
    u = a * p_prime / p
    d = a * (1.0 - p_prime) / (1.0 - p)

    df = np.exp(-r * dt)

    # Terminal option values.
    values = np.empty(num_steps + 1)

    s_low = stock_price * d**num_steps
    s = s_low

    ud_ratio = u / d

    for j in range(num_steps + 1):

        if is_call:
            values[j] = max(
                s - strike_price,
                0.0,
            )
        else:
            values[j] = max(
                strike_price - s,
                0.0,
            )

        s *= ud_ratio

    # Values required for Greeks.
    v10 = 0.0
    v11 = 0.0

    v20 = 0.0
    v21 = 0.0
    v22 = 0.0

    # Backward induction.
    for i in range(num_steps - 1, -1, -1):

        s_low /= d
        s = s_low

        for j in range(i + 1):

            hold = df * (
                p * values[j + 1]
                + (1.0 - p) * values[j]
            )

            if is_american:

                if is_call:
                    exercise = max(
                        s - strike_price,
                        0.0,
                    )
                else:
                    exercise = max(
                        strike_price - s,
                        0.0,
                    )

                values[j] = max(
                    exercise,
                    hold,
                )

            else:
                values[j] = hold

            s *= ud_ratio

        if i == 2:
            v20 = values[0]
            v21 = values[1]
            v22 = values[2]

        elif i == 1:
            v10 = values[0]
            v11 = values[1]

    price = values[0]

    # Greeks.
    s_dn = stock_price * d
    s_up = stock_price * u

    delta = (
        (v11 - v10)
        / (s_up - s_dn)
    )

    s_dd = stock_price * d * d
    s_ud = stock_price * d * u
    s_uu = stock_price * u * u

    delta_dn = (
        (v21 - v20)
        / (s_ud - s_dd)
    )

    delta_up = (
        (v22 - v21)
        / (s_uu - s_ud)
    )

    gamma = (
        (delta_up - delta_dn)
        / (s_up - s_dn)
    )

    theta = (
        (v21 - price)
        / (2.0 * dt)
    )

    return {
        "value": price,
        "delta": delta,
        "gamma": gamma,
        "theta": theta,
    }


########################################################################################
