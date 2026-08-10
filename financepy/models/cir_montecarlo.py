# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np
from numba import njit, float64, int64

from ..utils.helpers import label_to_string
from ..utils.global_types import CIRNumericalSchemeTypes


class CIRMonteCarlo:
    """Monte Carlo implementation of the Cox-Ingersoll-Ross model."""

    def __init__(self, a: float, b: float, sigma: float) -> None:
        _validate_model_params(a, b, sigma)
        self._a = a
        self._b = b
        self._sigma = sigma

    def __repr__(self) -> str:
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("Sigma", self._sigma)
        s += label_to_string("a", self._a)
        s += label_to_string("b", self._b)
        return s


@njit(cache=True)
def _validate_model_params(a: float, b: float, sigma: float) -> None:
    if a <= 0.0:
        raise ValueError("a must be positive.")
    if b < 0.0:
        raise ValueError("b must be non-negative.")
    if sigma < 0.0:
        raise ValueError("sigma must be non-negative.")


@njit(cache=True)
def _validate_time_params(t: float, dt: float) -> None:
    if t < 0.0:
        raise ValueError("t must be non-negative.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")


@njit(cache=True)
def _num_steps_and_dt(t: float, dt: float):
    if t == 0.0:
        return 0, dt

    num_steps = int(np.ceil(t / dt))
    dt_eff = t / num_steps
    return num_steps, dt_eff


@njit(cache=True)
def satisfies_feller(a: float, b: float, sigma: float) -> bool:
    _validate_model_params(a, b, sigma)
    return 2.0 * a * b >= sigma * sigma


@njit(cache=True)
def meanr(r0: float, a: float, b: float, t: float) -> float:
    _validate_model_params(a, b, 0.0)

    if r0 < 0.0:
        raise ValueError("r0 must be non-negative.")
    if t < 0.0:
        raise ValueError("t must be non-negative.")

    exp_at = np.exp(-a * t)
    return r0 * exp_at + b * (1.0 - exp_at)


@njit(cache=True)
def variancer(r0: float, a: float, b: float, sigma: float, t: float) -> float:
    _validate_model_params(a, b, sigma)

    if r0 < 0.0:
        raise ValueError("r0 must be non-negative.")
    if t < 0.0:
        raise ValueError("t must be non-negative.")

    exp_at = np.exp(-a * t)
    exp_2at = np.exp(-2.0 * a * t)

    vr = r0 * sigma * sigma * (exp_at - exp_2at) / a
    vr += b * sigma * sigma * ((1.0 - exp_at) ** 2) / (2.0 * a)
    return vr


@njit(
    float64(float64, float64, float64, float64, float64),
    cache=True,
)
def zero_price(r0: float, a: float, b: float, sigma: float, t: float) -> float:
    _validate_model_params(a, b, sigma)

    if r0 < 0.0:
        raise ValueError("r0 must be non-negative.")
    if t < 0.0:
        raise ValueError("t must be non-negative.")
    if t == 0.0:
        return 1.0

    if sigma == 0.0:
        integral_mean = b * t + (r0 - b) * (1.0 - np.exp(-a * t)) / a
        return np.exp(-integral_mean)

    sigma2 = sigma * sigma
    h = np.sqrt(a * a + 2.0 * sigma2)
    exp_ht = np.exp(h * t)

    denom = 2.0 * h + (a + h) * (exp_ht - 1.0)

    aa = (2.0 * h * np.exp((a + h) * t / 2.0) / denom) ** (2.0 * a * b / sigma2)

    bb = 2.0 * (exp_ht - 1.0) / denom

    return aa * np.exp(-r0 * bb)


@njit(
    float64(float64, float64, float64, float64, float64),
    cache=True,
)
def draw(rt: float, a: float, b: float, sigma: float, dt: float) -> float:
    _validate_model_params(a, b, sigma)
    _validate_time_params(dt, dt)

    if rt < 0.0:
        raise ValueError("rt must be non-negative.")

    if sigma == 0.0:
        return b + (rt - b) * np.exp(-a * dt)

    sigma2 = sigma * sigma
    exp_adt = np.exp(-a * dt)
    one_minus_exp = max(1.0 - exp_adt, 1.0e-16)

    d = 4.0 * a * b / sigma2
    ll = 4.0 * a * exp_adt * rt / (sigma2 * one_minus_exp)
    c = sigma2 * one_minus_exp / (4.0 * a)

    if d > 1.0:
        z = np.random.normal()
        x = np.random.chisquare(d - 1.0)
        r = c * (x + (z + np.sqrt(ll)) ** 2)
    else:
        n = np.random.poisson(ll / 2.0)
        x = np.random.chisquare(d + 2.0 * n)
        r = c * x

    return max(r, 0.0)


@njit(
    float64[:](float64, float64, float64, float64, float64, float64, int64, int64),
    cache=True,
)
def rate_path_mc(
    r0: float,
    a: float,
    b: float,
    sigma: float,
    t: float,
    dt: float,
    seed: int,
    scheme: int,
) -> np.ndarray:
    """Generate one CIR rate path from 0 to t inclusive."""

    _validate_model_params(a, b, sigma)
    _validate_time_params(t, dt)

    if r0 < 0.0:
        raise ValueError("r0 must be non-negative.")

    np.random.seed(seed)

    num_steps, dt_eff = _num_steps_and_dt(t, dt)

    rate_path = np.empty(num_steps + 1)
    rate_path[0] = r0

    if num_steps == 0:
        return rate_path

    if scheme == CIRNumericalSchemeTypes.EULER.value:

        sigmasqrt_dt = sigma * np.sqrt(dt_eff)
        r = r0

        for i_step in range(1, num_steps + 1):
            z = np.random.normal()
            r = r + a * (b - r) * dt_eff + z * sigmasqrt_dt * np.sqrt(max(r, 0.0))
            rate_path[i_step] = r

    elif scheme == CIRNumericalSchemeTypes.LOGNORMAL.value:

        x = np.exp(-a * dt_eff)
        y = 1.0 - x
        r = r0

        for i_step in range(1, num_steps + 1):
            mean = x * r + b * y
            var = sigma * sigma * y * (x * r + 0.5 * b * y) / a

            if sigma == 0.0 or var <= 0.0:
                r = mean
            else:
                sig = np.sqrt(np.log(1.0 + var / (mean * mean)))
                z = np.random.normal()
                r = mean * np.exp(-0.5 * sig * sig + sig * z)

            rate_path[i_step] = max(r, 0.0)

    elif scheme == CIRNumericalSchemeTypes.MILSTEIN.value:

        sigmasqrt_dt = sigma * np.sqrt(dt_eff)
        sigma2dt = sigma * sigma * dt_eff / 4.0
        r = r0

        for i_step in range(1, num_steps + 1):
            z = np.random.normal()
            r = r + a * (b - r) * dt_eff + z * sigmasqrt_dt * np.sqrt(max(r, 0.0))
            r = r + sigma2dt * (z * z - 1.0)
            rate_path[i_step] = r

    elif scheme == CIRNumericalSchemeTypes.KAHLJACKEL.value:

        if sigma == 0.0:
            r = r0
            for i_step in range(1, num_steps + 1):
                r = b + (r - b) * np.exp(-a * dt_eff)
                rate_path[i_step] = r

        else:
            bhat = b - sigma * sigma / (4.0 * a)
            sqrt_dt = np.sqrt(dt_eff)
            r = r0

            for i_step in range(1, num_steps + 1):
                z = np.random.normal()
                beta = z / sqrt_dt
                rootr = np.sqrt(max(r, 1.0e-12))
                c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt_eff / (4.0 * rootr)
                r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt_eff
                rate_path[i_step] = r

    elif scheme == CIRNumericalSchemeTypes.EXACT.value:

        r = r0

        for i_step in range(1, num_steps + 1):
            r = draw(r, a, b, sigma, dt_eff)
            rate_path[i_step] = r

    else:
        raise ValueError("Unknown CIR numerical scheme.")

    return rate_path


@njit(
    float64(
        float64,
        float64,
        float64,
        float64,
        float64,
        float64,
        int64,
        int64,
        int64,
    ),
    cache=True,
)
def zero_price_mc(
    r0: float,
    a: float,
    b: float,
    sigma: float,
    t: float,
    dt: float,
    num_paths: int,
    seed: int,
    scheme: int,
) -> float:
    """Estimate the CIR zero coupon price by Monte Carlo."""

    _validate_model_params(a, b, sigma)
    _validate_time_params(t, dt)

    if r0 < 0.0:
        raise ValueError("r0 must be non-negative.")

    if num_paths <= 0:
        raise ValueError("num_paths must be positive.")

    if t == 0.0:
        return 1.0

    np.random.seed(seed)

    num_steps, dt_eff = _num_steps_and_dt(t, dt)
    payoffs = np.empty(num_paths)

    for i_path in range(num_paths):

        r = r0
        rsum = r0

        if scheme == CIRNumericalSchemeTypes.EULER.value:

            sigmasqrt_dt = sigma * np.sqrt(dt_eff)

            for _ in range(1, num_steps + 1):
                r_prev = r
                z = np.random.normal()
                r = r + a * (b - r) * dt_eff + z * sigmasqrt_dt * np.sqrt(max(r, 0.0))
                rsum += r + r_prev

        elif scheme == CIRNumericalSchemeTypes.LOGNORMAL.value:

            x = np.exp(-a * dt_eff)
            y = 1.0 - x

            for _ in range(1, num_steps + 1):
                r_prev = r
                mean = x * r + b * y
                var = sigma * sigma * y * (x * r + 0.5 * b * y) / a

                if sigma == 0.0 or var <= 0.0:
                    r = mean
                else:
                    sig = np.sqrt(np.log(1.0 + var / (mean * mean)))
                    z = np.random.normal()
                    r = mean * np.exp(-0.5 * sig * sig + sig * z)

                r = max(r, 0.0)
                rsum += r + r_prev

        elif scheme == CIRNumericalSchemeTypes.MILSTEIN.value:

            sigmasqrt_dt = sigma * np.sqrt(dt_eff)
            sigma2dt = sigma * sigma * dt_eff / 4.0

            for _ in range(1, num_steps + 1):
                r_prev = r
                z = np.random.normal()
                r = r + a * (b - r) * dt_eff + z * sigmasqrt_dt * np.sqrt(max(r, 0.0))
                r = r + sigma2dt * (z * z - 1.0)
                rsum += r + r_prev

        elif scheme == CIRNumericalSchemeTypes.KAHLJACKEL.value:

            if sigma == 0.0:
                for _ in range(1, num_steps + 1):
                    r_prev = r
                    r = b + (r - b) * np.exp(-a * dt_eff)
                    rsum += r + r_prev

            else:
                bhat = b - sigma * sigma / (4.0 * a)
                sqrt_dt = np.sqrt(dt_eff)

                for _ in range(1, num_steps + 1):
                    r_prev = r
                    z = np.random.normal()
                    beta = z / sqrt_dt
                    rootr = np.sqrt(max(r, 1.0e-12))
                    c = 1.0 + (sigma * beta - 2.0 * a * rootr) * dt_eff / (4.0 * rootr)
                    r = r + (a * (bhat - r) + sigma * beta * rootr) * c * dt_eff
                    rsum += r + r_prev

        elif scheme == CIRNumericalSchemeTypes.EXACT.value:

            for _ in range(1, num_steps + 1):
                r_prev = r
                r = draw(r, a, b, sigma, dt_eff)
                rsum += r + r_prev

        else:
            raise ValueError("Unknown CIR numerical scheme.")

        payoffs[i_path] = np.exp(-0.5 * rsum * dt_eff)

    return np.mean(payoffs)
