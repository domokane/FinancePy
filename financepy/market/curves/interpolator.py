# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Refactored for clarity and faster vector evaluation.

from __future__ import annotations

from typing import Optional, Union

import numpy as np
from numba import njit, float64, int64
from scipy.interpolate import CubicSpline, PchipInterpolator

from ...utils.error import FinError
from ...utils.global_vars import G_SMALL
from ...utils.global_types import InterpTypes
from ...utils.tension_spline import TensionSpline



#    LINEAR_AVG_FWD_RATES = 11


# Groups used to make dispatch readable and cheap.
_FAST_NUMBA_TYPES = (
    InterpTypes.FLAT_FWD_RATES,
    InterpTypes.LINEAR_DISCOUNT,
    InterpTypes.LINEAR_ZERO_RATES,
)
_LOG_DF_SPLINE_TYPES = (
    InterpTypes.PCHIP_LOG_DISCOUNT,
    InterpTypes.NATCUBIC_LOG_DISCOUNT,
)
_ZERO_RATE_SPLINE_TYPES = (
    InterpTypes.PCHIP_ZERO_RATES,
    InterpTypes.FINCUBIC_ZERO_RATES,
    InterpTypes.NATCUBIC_ZERO_RATES,
    InterpTypes.TENSION_ZERO_RATES,
)

###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _zero_rates_from_dfs(times: np.ndarray, dfs: np.ndarray) -> np.ndarray:
    n = times.size
    out = np.empty(n, dtype=np.float64)
    for i in range(n):
        out[i] = -np.log(dfs[i]) / (times[i] + G_SMALL)
    if n > 1 and times[0] == 0.0:
        out[0] = out[1]
    return out


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _check_strictly_increasing(x: np.ndarray) -> bool:
    for i in range(1, x.size):
        if x[i] <= x[i - 1]:
            return False
    return True


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _find_interval(t: float, times: np.ndarray) -> int:
    """Return right index i such that t <= times[i].

    If t is beyond the last knot, return n. This mirrors the original code but
    uses binary search rather than a linear scan.
    """
    n = times.size
    lo = 0
    hi = n
    while lo < hi:
        mid = (lo + hi) // 2
        if times[mid] < t:
            lo = mid + 1
        else:
            hi = mid
    if lo == n:
        return n
    return lo


###############################################################################


@njit(
    float64(float64, float64[:], float64[:], int64),
    cache=True,
    fastmath=True,
    nogil=True,
)
def _uinterpolate(t: float, times: np.ndarray, dfs: np.ndarray, method: int) -> float:
    if t < 0.0:
        print(t, times, dfs, method)
        raise ValueError("Interpolation times must be non-negative.")
    if np.abs(t) < G_SMALL:
        return 1.0

    n = times.size
    if n == 0:
        return 1.0
    if n == 1:
        if np.abs(times[0]) < G_SMALL:
            return 1.0
        r = -np.log(dfs[0]) / times[0]
        return np.exp(-r * t)

    i = _find_interval(t, times)

    if method == InterpTypes.LINEAR_ZERO_RATES.value:
        if i == 1:
            r1 = -np.log(dfs[1]) / times[1]
            r2 = r1
            dt = times[1] - times[0]
            r = ((times[1] - t) * r1 + (t - times[0]) * r2) / dt
            return np.exp(-r * t)
        if i < n:
            r1 = -np.log(dfs[i - 1]) / times[i - 1]
            r2 = -np.log(dfs[i]) / times[i]
            dt = times[i] - times[i - 1]
            r = ((times[i] - t) * r1 + (t - times[i - 1]) * r2) / dt
            return np.exp(-r * t)
        r1 = -np.log(dfs[n - 1]) / times[n - 1]
        r2 = r1
        dt = times[n - 1] - times[n - 2]
        r = ((times[n - 1] - t) * r1 + (t - times[n - 2]) * r2) / dt
        return np.exp(-r * t)

    if method == InterpTypes.FLAT_FWD_RATES.value:
        if i == 0:
            return dfs[0]
        if i < n:
            left = i - 1
            right = i
        else:
            left = n - 2
            right = n - 1
        y1 = -np.log(dfs[left])
        y2 = -np.log(dfs[right])
        dt = times[right] - times[left]
        y = ((times[right] - t) * y1 + (t - times[left]) * y2) / dt
        return np.exp(-y)

    if method == InterpTypes.LINEAR_DISCOUNT.value:
        if i == 0:
            return dfs[0]

        if i < n:
            left = i - 1
            right = i
        else:
            left = n - 2
            right = n - 1

        dt = times[right] - times[left]

        df = (
            (times[right] - t) * dfs[left]
            + (t - times[left]) * dfs[right]
        ) / dt

        return df

    # Removed this method as it is not doing exactly what it claims
    # if method == InterpTypes.LINEAR_FWD_RATES.value:
    #     if i == 1:
    #         y = t * (-np.log(dfs[1] + 1.0e-10)) / (times[1] + 1.0e-10)
    #         return np.exp(-y)
    #     if i < n:
    #         f1 = -np.log(dfs[i - 1] / dfs[i - 2]) / (times[i - 1] - times[i - 2])
    #         f2 = -np.log(dfs[i] / dfs[i - 1]) / (times[i] - times[i - 1])
    #         dt = times[i] - times[i - 1]
    #         f = ((times[i] - t) * f1 + (t - times[i - 1]) * f2) / dt
    #         return dfs[i - 1] * np.exp(-f * (t - times[i - 1]))
    #     f = -np.log(dfs[n - 1] / dfs[n - 2]) / (times[n - 1] - times[n - 2])
    #     return dfs[n - 1] * np.exp(-f * (t - times[n - 1]))

    # if method == InterpTypes.LINEAR_AVG_FWD_RATES.value:
    #     if i == 1:
    #         y = t * (-np.log(dfs[1] + 1.0e-10)) / (times[1] + 1.0e-10)
    #         return np.exp(-y)

    #     if i < n:
    #         f1 = -np.log(dfs[i - 1] / dfs[i - 2]) / (times[i - 1] - times[i - 2])
    #         f2 = -np.log(dfs[i] / dfs[i - 1]) / (times[i] - times[i - 1])

    #         dt_total = times[i] - times[i - 1]
    #         dt = t - times[i - 1]

    #         slope = (f2 - f1) / dt_total
    #         integral = f1 * dt + 0.5 * slope * dt * dt

    #         return dfs[i - 1] * np.exp(-integral)

    #     f = -np.log(dfs[n - 1] / dfs[n - 2]) / (times[n - 1] - times[n - 2])
    #     return dfs[n - 1] * np.exp(-f * (t - times[n - 1]))

    raise ValueError("Invalid fast interpolation scheme.")


###############################################################################


@njit(
    float64[:](float64[:], float64[:], float64[:], int64),
    cache=True,
    fastmath=True,
    nogil=True,
)
def _vinterpolate(
    tvec: np.ndarray, times: np.ndarray, dfs: np.ndarray, method: int
) -> np.ndarray:
    out = np.empty(tvec.size, dtype=np.float64)
    for j in range(tvec.size):
        out[j] = _uinterpolate(tvec[j], times, dfs, method)
    return out


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _build_linear_onf_curve(times: np.ndarray, dfs: np.ndarray):
    tmp_t = np.empty(times.size + 1, dtype=np.float64)
    tmp_r = np.empty(times.size + 1, dtype=np.float64)
    count = 0
    prev_df = 1.0

    for j in range(times.size):
        t = times[j]
        df = dfs[j]
        if t == 0.0:
            prev_df = df
            continue

        if count == 0:
            r = -np.log(df) / t
            tmp_t[0] = 0.0
            tmp_r[0] = r
            tmp_t[1] = t
            tmp_r[1] = r
            count = 2
        else:
            prev_t = tmp_t[count - 1]
            prev_r = tmp_r[count - 1]
            onf_int = -np.log(df / prev_df)
            r = 2.0 * onf_int / (t - prev_t) - prev_r
            tmp_t[count] = t
            tmp_r[count] = r
            count += 1
        prev_df = df

    if count == 0:
        return np.array([0.0, 0.1]), np.array([0.0, 0.0])
    return tmp_t[:count], tmp_r[:count]


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _linear_onf_cumulative_integrals(
    onf_times: np.ndarray, onf_rates: np.ndarray
) -> np.ndarray:
    n = onf_times.size
    cum = np.zeros(n, dtype=np.float64)
    for i in range(1, n):
        dt = onf_times[i] - onf_times[i - 1]
        cum[i] = cum[i - 1] + 0.5 * dt * (onf_rates[i - 1] + onf_rates[i])
    return cum


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _linear_onf_integral(
    t: float, onf_times: np.ndarray, onf_rates: np.ndarray, cum: np.ndarray
) -> float:
    if t <= 0.0:
        return 0.0
    n = onf_times.size
    i = _find_interval(t, onf_times)

    if i == 0:
        return t * onf_rates[0]
    if i >= n:
        return cum[n - 1] + (t - onf_times[n - 1]) * onf_rates[n - 1]

    t0 = onf_times[i - 1]
    t1 = onf_times[i]
    r0 = onf_rates[i - 1]
    r1 = onf_rates[i]
    dt = t - t0
    slope = (r1 - r0) / (t1 - t0)
    return cum[i - 1] + r0 * dt + 0.5 * slope * dt * dt


###############################################################################


@njit(cache=True, fastmath=True, nogil=True)
def _linear_onf_dfs(
    tvec: np.ndarray,
    onf_times: np.ndarray,
    onf_rates: np.ndarray,
    cum: np.ndarray,
) -> np.ndarray:
    out = np.empty(tvec.size, dtype=np.float64)
    for i in range(tvec.size):
        if tvec[i] < 0.0:
            raise ValueError("Interpolation times must be non-negative.")
        out[i] = np.exp(-_linear_onf_integral(tvec[i], onf_times, onf_rates, cum))
    return out


###############################################################################


def interpolate(
    t: Union[float, np.ndarray],
    times: np.ndarray,
    dfs: np.ndarray,
    method: int,
):
    """Fast stand-alone discount-factor interpolation.

    This compatibility wrapper supports the original three Numba-backed methods:
    flat forwards, linear forwards and linear zero rates.
    """
    times = np.asarray(times, dtype=np.float64)
    dfs = np.asarray(dfs, dtype=np.float64)
    if isinstance(t, (float, np.float64)):
        return _uinterpolate(float(t), times, dfs, int(method))
    if isinstance(t, np.ndarray):
        return _vinterpolate(np.asarray(t, dtype=np.float64), times, dfs, int(method))
    raise FinError(f"Unknown input type {type(t)}")


###############################################################################


class Interpolator:
    """Discount-factor interpolator.

    The class separates fitting from evaluation:

    * ``fit`` validates and stores the market knots, and builds any expensive
      spline representation once.
    * ``interpolate`` performs the cheapest possible evaluation path for the
      selected interpolation type.
    """

    def __init__(self, interpolator_type: InterpTypes, **kwargs: dict):
        self._interp_type = interpolator_type
        self._method = interpolator_type.value
        self._times: Optional[np.ndarray] = None
        self._dfs: Optional[np.ndarray] = None
        self._interp_fn = None
        self._onf_times: Optional[np.ndarray] = None
        self._onf_rates: Optional[np.ndarray] = None
        self._onf_integrals: Optional[np.ndarray] = None
        self._optional_interp_params = kwargs

    def fit(self, times: np.ndarray, dfs: np.ndarray) -> None:
        times = np.asarray(times, dtype=np.float64)
        dfs = np.asarray(dfs, dtype=np.float64)

        if times.ndim != 1 or dfs.ndim != 1:
            raise FinError("times and dfs must be one-dimensional arrays.")
        if times.size != dfs.size:
            raise FinError("times and dfs must have the same length.")
        if np.any(times < 0.0):
            raise FinError("Curve times must be non-negative.")
        if np.any(dfs <= 0.0):
            raise FinError("Discount factors must be positive.")
        if times.size > 1 and not _check_strictly_increasing(times):
            raise FinError("times must be strictly increasing.")

        self._times = times
        self._dfs = dfs
        self._interp_fn = None
        self._onf_times = None
        self._onf_rates = None
        self._onf_integrals = None

        if times.size <= 1 or self._interp_type in _FAST_NUMBA_TYPES:
            return

        if self._interp_type == InterpTypes.PCHIP_LOG_DISCOUNT:
            self._interp_fn = PchipInterpolator(times, np.log(dfs), extrapolate=True)
            return

        if self._interp_type == InterpTypes.PCHIP_ZERO_RATES:
            self._interp_fn = PchipInterpolator(
                times, _zero_rates_from_dfs(times, dfs), extrapolate=True
            )
            return

        if self._interp_type == InterpTypes.FINCUBIC_ZERO_RATES:
            self._interp_fn = CubicSpline(
                times,
                _zero_rates_from_dfs(times, dfs),
                bc_type=((2, 0.0), (1, 0.0)),
                extrapolate=True,
            )
            return

        if self._interp_type == InterpTypes.NATCUBIC_LOG_DISCOUNT:
            self._interp_fn = CubicSpline(
                times, np.log(dfs), bc_type="natural", extrapolate=True
            )
            return

        if self._interp_type == InterpTypes.NATCUBIC_ZERO_RATES:
            self._interp_fn = CubicSpline(
                times,
                _zero_rates_from_dfs(times, dfs),
                bc_type="natural",
                extrapolate=True,
            )
            return

        if self._interp_type == InterpTypes.LINEAR_ONFWD_RATES:
            self._onf_times, self._onf_rates = _build_linear_onf_curve(times, dfs)
            self._onf_integrals = _linear_onf_cumulative_integrals(
                self._onf_times, self._onf_rates
            )
            return

        if self._interp_type == InterpTypes.TENSION_ZERO_RATES:
            sigma = self._optional_interp_params.get("sigma", 1.0)
            self._interp_fn = TensionSpline(
                times, _zero_rates_from_dfs(times, dfs), sigma=sigma
            )
            return

        raise FinError(f"Unknown interpolation type {self._interp_type}")

    def interpolate(self, t: Union[float, np.ndarray]):
        if self._times is None or self._dfs is None:
            raise FinError("Interpolator has not been fitted.")

        scalar = isinstance(t, (float, np.float64, int, np.integer))
        if scalar:
            tvec = np.array([float(t)], dtype=np.float64)
        elif isinstance(t, np.ndarray):
            tvec = np.asarray(t, dtype=np.float64)
        else:
            raise FinError(f"t is not a recognized type: {type(t)}")

        if np.any(tvec < 0.0):
            raise FinError("Interpolation times must be non-negative.")

        if self._times.size <= 1:
            out = _vinterpolate(tvec, self._times, self._dfs, self._method)
        elif self._interp_type in _FAST_NUMBA_TYPES:
            out = _vinterpolate(tvec, self._times, self._dfs, self._method)
        elif self._interp_type in _LOG_DF_SPLINE_TYPES:
            out = np.exp(self._interp_fn(tvec))
        elif self._interp_type in _ZERO_RATE_SPLINE_TYPES:
            out = np.exp(-tvec * self._interp_fn(tvec))
        elif self._interp_type == InterpTypes.LINEAR_ONFWD_RATES:
            out = _linear_onf_dfs(
                tvec, self._onf_times, self._onf_rates, self._onf_integrals
            )
        else:
            raise FinError(f"Unknown interpolation type {self._interp_type}")

        if scalar:
            return float(out[0])
        else:
            return out

    @classmethod
    def suitable_for_bootstrap(cls, interp_type: InterpTypes) -> bool:
        return interp_type in {
            InterpTypes.FLAT_FWD_RATES,
            InterpTypes.LINEAR_DISCOUNT,
            InterpTypes.LINEAR_ZERO_RATES,
            InterpTypes.LINEAR_ONFWD_RATES,
        }
