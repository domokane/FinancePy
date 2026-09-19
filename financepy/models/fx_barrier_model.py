from math import erf, exp, log, sqrt

import numba as nb
import numpy as np

_SQRT_2 = sqrt(2.0)


@nb.njit(cache=True)
def _norm_cdf(x: float) -> float:
    """Numba-compatible scalar normal CDF."""
    return 0.5 * (1.0 + erf(x / _SQRT_2))


@nb.njit(cache=True)
def _fx_barrier_price_scalar(
    s0: float,
    k: float,
    h: float,
    t: float,
    df: float,
    dq: float,
    volatility: float,
    num_obs_per_year: int,
    barrier_type: int,
) -> float:
    """Price one FX barrier option for one spot value."""

    # Preserve the existing implementation's minimum volatility,
    # but apply it before sigma_root_t and v2 are calculated.
    sigma = max(abs(volatility), 1.0e-5)

    sqrt_t = sqrt(t)
    sigma_root_t = sigma * sqrt_t
    v2 = sigma * sigma

    r_d = -log(df) / t
    r_f = -log(dq) / t
    mu = r_d - r_f

    ln_s0_k = log(s0 / k)

    d1 = (ln_s0_k + (mu + 0.5 * v2) * t) / sigma_root_t
    d2 = (ln_s0_k + (mu - 0.5 * v2) * t) / sigma_root_t

    call_value = (
        s0 * dq * _norm_cdf(d1)
        - k * df * _norm_cdf(d2)
    )

    put_value = (
        k * df * _norm_cdf(-d2)
        - s0 * dq * _norm_cdf(-d1)
    )

    # Enum values currently used by FXBarrierTypes:
    # 1 DOWN_AND_OUT_CALL
    # 2 DOWN_AND_IN_CALL
    # 3 UP_AND_OUT_CALL
    # 4 UP_AND_IN_CALL
    # 5 UP_AND_OUT_PUT
    # 6 UP_AND_IN_PUT
    # 7 DOWN_AND_OUT_PUT
    # 8 DOWN_AND_IN_PUT

    if barrier_type == 1 and s0 <= h:
        return 0.0
    if barrier_type == 3 and s0 >= h:
        return 0.0
    if barrier_type == 5 and s0 >= h:
        return 0.0
    if barrier_type == 7 and s0 <= h:
        return 0.0

    if barrier_type == 2 and s0 <= h:
        return call_value
    if barrier_type == 4 and s0 >= h:
        return call_value
    if barrier_type == 6 and s0 >= h:
        return put_value
    if barrier_type == 8 and s0 <= h:
        return put_value

    # The existing expression:
    #
    # num_observations = t * num_obs_per_year
    # sqrt(t / num_observations)
    #
    # simplifies to sqrt(1 / num_obs_per_year).
    monitoring_adjustment = (
        0.5826 * sigma / sqrt(float(num_obs_per_year))
    )

    if barrier_type in (1, 2, 7, 8):
        h *= exp(-monitoring_adjustment)
    elif barrier_type in (3, 4, 5, 6):
        h *= exp(monitoring_adjustment)
    else:
        raise ValueError("Unknown barrier option type")

    ll = (mu + 0.5 * v2) / v2

    y = (
        log(h * h / (s0 * k)) / sigma_root_t
        + ll * sigma_root_t
    )

    x1 = (
        log(s0 / h) / sigma_root_t
        + ll * sigma_root_t
    )

    y1 = (
        log(h / s0) / sigma_root_t
        + ll * sigma_root_t
    )

    h_over_s = h / s0
    power_1 = h_over_s ** (2.0 * ll)
    power_2 = h_over_s ** (2.0 * ll - 2.0)

    if barrier_type == 1:  # Down-and-out call
        if h >= k:
            return (
                s0 * dq * _norm_cdf(x1)
                - k * df * _norm_cdf(x1 - sigma_root_t)
                - s0 * dq * power_1 * _norm_cdf(y1)
                + k * df * power_2 * _norm_cdf(y1 - sigma_root_t)
            )

        call_in = (
            s0 * dq * power_1 * _norm_cdf(y)
            - k * df * power_2 * _norm_cdf(y - sigma_root_t)
        )
        return call_value - call_in

    if barrier_type == 2:  # Down-and-in call
        if h <= k:
            return (
                s0 * dq * power_1 * _norm_cdf(y)
                - k * df * power_2 * _norm_cdf(y - sigma_root_t)
            )

        call_out = (
            s0 * dq * _norm_cdf(x1)
            - k * df * _norm_cdf(x1 - sigma_root_t)
            - s0 * dq * power_1 * _norm_cdf(y1)
            + k * df * power_2 * _norm_cdf(y1 - sigma_root_t)
        )
        return call_value - call_out

    if barrier_type == 4:  # Up-and-in call
        if h >= k:
            return (
                s0 * dq * _norm_cdf(x1)
                - k * df * _norm_cdf(x1 - sigma_root_t)
                - s0 * dq * power_1
                * (_norm_cdf(-y) - _norm_cdf(-y1))
                + k * df * power_2
                * (
                    _norm_cdf(-y + sigma_root_t)
                    - _norm_cdf(-y1 + sigma_root_t)
                )
            )

        return call_value

    if barrier_type == 3:  # Up-and-out call
        if h <= k:
            return 0.0

        call_in = (
            s0 * dq * _norm_cdf(x1)
            - k * df * _norm_cdf(x1 - sigma_root_t)
            - s0 * dq * power_1
            * (_norm_cdf(-y) - _norm_cdf(-y1))
            + k * df * power_2
            * (
                _norm_cdf(-y + sigma_root_t)
                - _norm_cdf(-y1 + sigma_root_t)
            )
        )
        return call_value - call_in

    if barrier_type == 6:  # Up-and-in put
        if h > k:
            return (
                -s0 * dq * power_1 * _norm_cdf(-y)
                + k * df * power_2
                * _norm_cdf(-y + sigma_root_t)
            )

        put_out = (
            -s0 * dq * _norm_cdf(-x1)
            + k * df * _norm_cdf(-x1 + sigma_root_t)
            + s0 * dq * power_1 * _norm_cdf(-y1)
            - k * df * power_2
            * _norm_cdf(-y1 + sigma_root_t)
        )
        return put_value - put_out

    if barrier_type == 5:  # Up-and-out put
        if h >= k:
            put_in = (
                -s0 * dq * power_1 * _norm_cdf(-y)
                + k * df * power_2
                * _norm_cdf(-y + sigma_root_t)
            )
            return put_value - put_in

        return (
            -s0 * dq * _norm_cdf(-x1)
            + k * df * _norm_cdf(-x1 + sigma_root_t)
            + s0 * dq * power_1 * _norm_cdf(-y1)
            - k * df * power_2
            * _norm_cdf(-y1 + sigma_root_t)
        )

    if barrier_type == 7:  # Down-and-out put
        if h >= k:
            return 0.0

        put_in = (
            -s0 * dq * _norm_cdf(-x1)
            + k * df * _norm_cdf(-x1 + sigma_root_t)
            + s0 * dq * power_1
            * (_norm_cdf(y) - _norm_cdf(y1))
            - k * df * power_2
            * (
                _norm_cdf(y - sigma_root_t)
                - _norm_cdf(y1 - sigma_root_t)
            )
        )
        return put_value - put_in

    if barrier_type == 8:  # Down-and-in put
        if h >= k:
            return put_value

        return (
            -s0 * dq * _norm_cdf(-x1)
            + k * df * _norm_cdf(-x1 + sigma_root_t)
            + s0 * dq * power_1
            * (_norm_cdf(y) - _norm_cdf(y1))
            - k * df * power_2
            * (
                _norm_cdf(y - sigma_root_t)
                - _norm_cdf(y1 - sigma_root_t)
            )
        )

    raise ValueError("Unknown barrier option type")


@nb.njit(cache=True, parallel=True)
def _fx_barrier_price_array(
    spots: np.ndarray,
    k: float,
    h: float,
    t: float,
    df: float,
    dq: float,
    volatility: float,
    num_obs_per_year: int,
    barrier_type: int,
) -> np.ndarray:
    """Price a flattened array of spot values."""

    values = np.empty(spots.size, dtype=np.float64)

    for index in nb.prange(spots.size):
        values[index] = _fx_barrier_price_scalar(
            spots[index],
            k,
            h,
            t,
            df,
            dq,
            volatility,
            num_obs_per_year,
            barrier_type,
        )

    return values


def fx_barrier_value(
    spot_fx_rate,
    strike_fx_rate: float,
    barrier_level: float,
    time_to_expiry: float,
    domestic_df: float,
    foreign_df: float,
    volatility: float,
    num_obs_per_year: int,
    barrier_type: int,
):
    """Price scalar or array-valued FX spots while preserving input shape."""

    spots = np.asarray(spot_fx_rate, dtype=np.float64)
    scalar_input = spots.ndim == 0

    if np.any(spots <= 0.0):
        raise ValueError("Spot FX rates must be positive")

    if strike_fx_rate <= 0.0:
        raise ValueError("Strike FX rate must be positive")

    if barrier_level <= 0.0:
        raise ValueError("Barrier level must be positive")

    if time_to_expiry <= 0.0:
        raise ValueError("Time to expiry must be positive")

    if num_obs_per_year <= 0:
        raise ValueError("num_obs_per_year must be positive")

    flat_spots = np.ascontiguousarray(spots.reshape(-1))

    flat_values = _fx_barrier_price_array(
        flat_spots,
        float(strike_fx_rate),
        float(barrier_level),
        float(time_to_expiry),
        float(domestic_df),
        float(foreign_df),
        float(volatility),
        int(num_obs_per_year),
        int(barrier_type),
    )

    if scalar_input:
        return float(flat_values[0])

    return flat_values.reshape(spots.shape)
