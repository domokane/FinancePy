# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from __future__ import annotations

from enum import Enum
from typing import Optional, Union

import numpy as np
from numba import njit, float64, int64

from scipy.interpolate import PchipInterpolator
from scipy.interpolate import CubicSpline
from scipy.interpolate import InterpolatedUnivariateSpline

from ...utils.error import FinError
from ...utils.global_vars import G_SMALL
from ...utils.tension_spline import TensionSpline

########################################################################################


class InterpTypes(Enum):
    FLAT_FWD_RATES = 1
    LINEAR_FWD_RATES = 2
    LINEAR_ZERO_RATES = 4
    FINCUBIC_ZERO_RATES = 7
    NATCUBIC_LOG_DISCOUNT = 8
    NATCUBIC_ZERO_RATES = 9
    PCHIP_ZERO_RATES = 10
    PCHIP_LOG_DISCOUNT = 11
    LINEAR_ONFWD_RATES = 21
    TENSION_ZERO_RATES = 22


# TODO: GET RID OF THIS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########################################################################################


def interpolate(
    t: Union[float, np.ndarray],  # time or array of times
    times: np.ndarray,  # Vector of times on grid
    dfs: np.ndarray,  # Vector of discount factors
    method: int,
):  # Interpolation method which is value of enum
    """Fast interpolation of discount factors at time x given discount factors
    at times provided using one of the methods in the enum InterpTypes. The
    value of x can be an array so that the function is vectorised."""

    if isinstance(t, (float, np.float64)):

        if t < 0.0:
            print(t)
            raise FinError("Interpolate times must all be >= 0")

        u = _uinterpolate(t, times, dfs, method)
        return u

    elif isinstance(t, np.ndarray):

        if np.any(t < 0.0):
            print(t)
            raise FinError("Interpolate times must all be >= 0")

        v = _vinterpolate(t, times, dfs, method)

        return v

    raise FinError("Unknown input type" + type(t))


@njit(
    float64(float64, float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
    nogil=True,
)

########################################################################################


def _uinterpolate(t, times, dfs, method):
    """Return the interpolated value of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate."""

    small = 1e-10
    num_points = times.size

    if t == times[0]:
        return dfs[0]

    i = 0
    while times[i] < t and i < num_points - 1:
        i = i + 1

    if t > times[i]:
        i = num_points

    yvalue = 0.0

    # linear interpolation of y(x)

    if method == InterpTypes.LINEAR_ZERO_RATES.value:

        if i == 1:
            r1 = -np.log(dfs[i]) / times[i]
            r2 = r1
            dt = times[i] - times[i - 1]
            rvalue = ((times[i] - t) * r1 + (t - times[i - 1]) * r2) / dt
            yvalue = np.exp(-rvalue * t)
        elif i < num_points:
            r1 = -np.log(dfs[i - 1]) / times[i - 1]
            r2 = -np.log(dfs[i]) / times[i]
            dt = times[i] - times[i - 1]
            rvalue = ((times[i] - t) * r1 + (t - times[i - 1]) * r2) / dt
            yvalue = np.exp(-rvalue * t)
        else:
            r1 = -np.log(dfs[i - 1]) / times[i - 1]
            r2 = -np.log(dfs[i - 1]) / times[i - 1]
            dt = times[i - 1] - times[i - 2]
            rvalue = ((times[i - 1] - t) * r1 + (t - times[i - 2]) * r2) / dt
            yvalue = np.exp(-rvalue * t)

        return yvalue

    # linear interpolation of log(y(x)) which means the linear interpolation of
    # continuously compounded zero rates in the case of discount discount
    # This is also FLAT FORWARDS

    elif method == InterpTypes.FLAT_FWD_RATES.value:

        if i == 1:
            rt1 = -np.log(dfs[i - 1])
            rt2 = -np.log(dfs[i])
            dt = times[i] - times[i - 1]
            rtvalue = ((times[i] - t) * rt1 + (t - times[i - 1]) * rt2) / dt
            yvalue = np.exp(-rtvalue)
        elif i < num_points:
            rt1 = -np.log(dfs[i - 1])
            rt2 = -np.log(dfs[i])
            dt = times[i] - times[i - 1]
            rtvalue = ((times[i] - t) * rt1 + (t - times[i - 1]) * rt2) / dt
            yvalue = np.exp(-rtvalue)
        else:
            rt1 = -np.log(dfs[i - 2])
            rt2 = -np.log(dfs[i - 1])
            dt = times[i - 1] - times[i - 2]
            rtvalue = (
                (times[i - 1] - t) * rt1 + (t - times[i - 2]) * rt2
            ) / dt
            yvalue = np.exp(-rtvalue)

        return yvalue

    elif method == InterpTypes.LINEAR_FWD_RATES.value:

        if i == 1:
            y2 = -np.log(dfs[i] + small)
            yvalue = t * y2 / (times[i] + small)
            yvalue = np.exp(-yvalue)
        elif i < num_points:
            # If you get a math domain error it is because you need negativ
            fwd1 = -np.log(dfs[i - 1] / dfs[i - 2]) / (
                times[i - 1] - times[i - 2]
            )
            fwd2 = -np.log(dfs[i] / dfs[i - 1]) / (times[i] - times[i - 1])
            dt = times[i] - times[i - 1]
            fwd = ((times[i] - t) * fwd1 + (t - times[i - 1]) * fwd2) / dt
            yvalue = dfs[i - 1] * np.exp(-fwd * (t - times[i - 1]))
        else:
            fwd = -np.log(dfs[i - 1] / dfs[i - 2]) / (
                times[i - 1] - times[i - 2]
            )
            yvalue = dfs[i - 1] * np.exp(-fwd * (t - times[i - 1]))

        return yvalue

    else:
        print(method)
        raise FinError("Invalid interpolation scheme.")


@njit(
    float64[:](float64[:], float64[:], float64[:], int64),
    fastmath=True,
    cache=True,
    nogil=True,
)

########################################################################################


def _vinterpolate(x_values, x_vector, dfs, method):
    """Return the interpolated values of y given x and a vector of x and y.
    The values of x must be monotonic and increasing. The different schemes for
    interpolation are linear in y (as a function of x), linear in log(y) and
    piecewise flat in the continuously compounded forward y rate."""

    n = x_values.size
    yvalues = np.empty(n)
    for i in range(0, n):
        yvalues[i] = _uinterpolate(x_values[i], x_vector, dfs, method)

    return yvalues


########################################################################################


class Interpolator:
    """Interpolator class for curve fitting and value interpolation."""

    ####################################################################################

    def __init__(self, interpolator_type: InterpTypes, **kwargs: dict):

        self._interp_type = interpolator_type
        self._interp_fn = None
        self._times = None
        self._dfs = None
        self._refit_curve = False
        self._optional_interp_params = kwargs

    ####################################################################################

    def fit(self, times: np.ndarray, dfs: np.ndarray):
        """
        Fit the interpolator to the provided time and discount factor arrays.
        """
        self._times = np.asarray(times, dtype=np.float64)
        self._dfs = np.asarray(dfs, dtype=np.float64)

        if len(self._times) <= 1:
            self._interp_fn = None
            return

        if self._interp_type == InterpTypes.PCHIP_LOG_DISCOUNT:

            log_dfs = np.log(self._dfs)
            self._interp_fn = PchipInterpolator(self._times, log_dfs)

        elif self._interp_type == InterpTypes.PCHIP_ZERO_RATES:

            g_small_vector = np.ones(len(self._times)) * G_SMALL
            zero_rates = -np.log(self._dfs) / (self._times + g_small_vector)

            if self._times[0] == 0.0:
                zero_rates[0] = zero_rates[1]

            self._interp_fn = PchipInterpolator(self._times, zero_rates)

        # if self._interp_type == InterpTypes.FINCUBIC_LOG_DISCOUNT:

        #     """ Second derivatives at left is zero and first derivative at
        #     right is clamped to zero. """
        #     log_dfs = np.log(self.dfs)
        #     self._interp_fn = CubicSpline(self.times, log_dfs,
        #                                  bc_type=((2, 0.0), (1, 0.0)))

        elif self._interp_type == InterpTypes.FINCUBIC_ZERO_RATES:

            # Second derivatives at left is zero and first derivative at
            # right is clamped to zero.
            g_small_vector = np.ones(len(self._times)) * G_SMALL
            zero_rates = -np.log(self._dfs) / (self._times + g_small_vector)

            if self._times[0] == 0.0:
                zero_rates[0] = zero_rates[1]

            self._interp_fn = CubicSpline(
                self._times, zero_rates, bc_type=((2, 0.0), (1, 0.0))
            )

        elif self._interp_type == InterpTypes.NATCUBIC_LOG_DISCOUNT:

            # Second derivatives are clamped to zero at end points
            log_dfs = np.log(self._dfs)
            self._interp_fn = CubicSpline(
                self._times, log_dfs, bc_type="natural"
            )

        elif self._interp_type == InterpTypes.NATCUBIC_ZERO_RATES:

            # Second derivatives are clamped to zero at end points
            g_small_vector = np.ones(len(self._times)) * G_SMALL
            zero_rates = -np.log(self._dfs) / (self._times + g_small_vector)

            if self._times[0] == 0.0:
                zero_rates[0] = zero_rates[1]

            self._interp_fn = CubicSpline(
                self._times, zero_rates, bc_type="natural"
            )

        #        elif self._interp_type  == InterpTypes.LINEAR_LOG_DISCOUNT:
        #            log_dfs = np.log(self.dfs)
        #            self._interp_fn = interp1d(self.times, log_dfs,
        #                                      fill_value="extrapolate")

        elif self._interp_type == InterpTypes.LINEAR_ONFWD_RATES:
            onf_times = []
            onf_rates = []
            prev_df = 1.0
            for t, df in zip(self._times, self._dfs):
                if t == 0.0:
                    continue
                if len(onf_times) == 0:
                    onfr = -np.log(df) / t
                    onf_times = [0.0, t]
                    onf_rates = [onfr, onfr]
                else:
                    prev_t = onf_times[-1]
                    prev_r = onf_rates[-1]
                    fwd_df = df / prev_df
                    onf_int = -np.log(fwd_df)
                    r = 2 * onf_int / (t - prev_t) - prev_r

                    onf_times.append(t)
                    onf_rates.append(r)

                prev_df = df

            if len(onf_times) == 0:
                self._interp_fn = InterpolatedUnivariateSpline(
                    [0.0, 0.1], [0.0, 0.0], k=1, ext=3
                )
            else:
                self._interp_fn = InterpolatedUnivariateSpline(
                    onf_times, onf_rates, k=1, ext=3
                )

        elif self._interp_type == InterpTypes.TENSION_ZERO_RATES:
            tension_sigma = self._optional_interp_params.get("sigma", 1.0)
            g_small_vector = np.ones(len(self._times)) * G_SMALL
            zero_rates = -np.log(self._dfs) / (self._times + g_small_vector)

            if self._times[0] == 0.0:
                zero_rates[0] = zero_rates[1]

            self._interp_fn = TensionSpline(
                self._times, zero_rates, sigma=tension_sigma
            )

    ####################################################################################

    def interpolate(self, t: float):
        """Interpolation of discount factors at time x given discount factors
        at times provided using one of the methods in the enum InterpTypes.
        The value of x can be an array so that the function is vectorised."""

        if self._dfs is None:
            raise FinError("Dfs have not been set.")

        if isinstance(t, (float, np.float64)):

            if t < 0.0:
                raise FinError("Interpolate times must all be >= 0")

            if np.abs(t) < G_SMALL:
                return 1.0

            scalar = True
            tvec = np.array([t])

        elif isinstance(t, np.ndarray):

            if np.any(t < 0.0):
                print(t)
                raise FinError("Interpolate times must all be >= 0")
            scalar = False
            tvec = t

        else:
            raise FinError("t is not a recognized type")

        # SPECIAL CASES: len(times) <= 1
        if self._interp_fn is None:
            n = len(self._times)

            if n == 0:
                # Empty fit → flat DF = 1.0
                out = np.full_like(tvec, 1.0)

            elif n == 1:
                t0 = self._times[0]
                df0 = self._dfs[0]

                if abs(t0) < G_SMALL:
                    # Single point at origin
                    out = np.full_like(tvec, 1.0)
                else:
                    # Single point NOT at origin → constant zero rate from t=0
                    rate = -np.log(df0) / t0
                    out = np.exp(-rate * tvec)

            else:
                # Should not happen
                out = _vinterpolate(
                    tvec, self._times, self._dfs, self._interp_type.value
                )

            return out[0] if scalar else out

        # NORMAL FITTED CASES

        # if self._interp_fn is None:
        #     # Single point at origin (common during bootstrap edge cases)
        #     if len(self._times) == 1 and abs(self._times[0]) < G_SMALL:
        #         if isinstance(t, (float, np.float64)):
        #             return 1.0
        #         else:
        #             return np.full_like(t, 1.0)

        #     # Otherwise fall back to the fast numba interpolator
        #     out = _vinterpolate(
        #         tvec, self._times, self._dfs, self._interp_type.value
        #     )

        if self._interp_type in [
            InterpTypes.PCHIP_LOG_DISCOUNT,
            InterpTypes.NATCUBIC_LOG_DISCOUNT,
        ]:

            out = np.exp(self._interp_fn(tvec))

        elif self._interp_type in [
            InterpTypes.PCHIP_ZERO_RATES,
            InterpTypes.FINCUBIC_ZERO_RATES,
            InterpTypes.NATCUBIC_ZERO_RATES,
            InterpTypes.TENSION_ZERO_RATES,
        ]:

            out = np.exp(-tvec * self._interp_fn(tvec))

        elif self._interp_type == InterpTypes.LINEAR_ONFWD_RATES:

            def true_integral(spline, t):

                last_t = spline.get_knots()[-1]
                i1 = spline.integral(0.0, min(t, last_t))
                i2 = (t - min(t, last_t)) * spline(last_t)
                return i1 + i2

            log_dfs = [-true_integral(self._interp_fn, t) for t in tvec]
            out = np.exp(log_dfs)

        else:

            out = _vinterpolate(
                tvec, self._times, self._dfs, self._interp_type.value
            )

        if scalar:
            return out[0]

        return out

    @classmethod

    ####################################################################################

    def suitable_for_bootstrap(cls, interp_type):
        """Check if the interpolation type is suitable for bootstrapping."""
        is_suitable = {
            InterpTypes.FLAT_FWD_RATES: True,
            InterpTypes.LINEAR_FWD_RATES: True,
            InterpTypes.LINEAR_ZERO_RATES: True,
            InterpTypes.FINCUBIC_ZERO_RATES: False,
            InterpTypes.NATCUBIC_LOG_DISCOUNT: False,
            InterpTypes.NATCUBIC_ZERO_RATES: False,
            InterpTypes.PCHIP_ZERO_RATES: False,
            InterpTypes.PCHIP_LOG_DISCOUNT: False,
            InterpTypes.LINEAR_ONFWD_RATES: True,
            InterpTypes.TENSION_ZERO_RATES: False,
        }

        return is_suitable[interp_type]
